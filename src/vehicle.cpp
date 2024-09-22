#include "constants.hpp"
#include "gnc.hpp"
#include "vehicle.hpp"

#include <exception>
#include <iostream>
#include <stdio.h>

InertialProperties::InertialProperties(double mass, std::array<double, 4> MOI, Eigen::Vector3d COM) :
    _mass_empty(mass), _MOI_empty(std::move(MOI)), _COM_empty(std::move(COM)) 
{
    _MOI = MOI;
    _COM = COM;
    _MOI_delta.fill(0.0);
    _COM_delta.setZero();
}

InertialProperties::InertialProperties(double mass_empty, std::array<double, 4> MOI_empty, Eigen::Vector3d COM_empty,
    double mass_full, const std::array<double, 4>& MOI_full, const Eigen::Vector3d& COM_full) : 
        _mass_empty(mass_empty), _MOI_empty(std::move(MOI_empty)), _COM_empty(std::move(COM_empty)) 
{
    if(mass_full < mass_empty)
    {
        throw new std::runtime_error("mass full must be larger than mass empty");
    }
    double dm = mass_full - mass_empty;
    for(auto i = 0u; i < 3; i++)
    {
        _MOI_delta[i] = (MOI_full[i] - MOI_empty[i])/dm;
        _COM_delta(i) = (COM_full[i] - COM_empty[i])/dm;
    }
    _MOI_delta[3] = (MOI_full[3] - MOI_empty[3])/dm;
}

void InertialProperties::set_mass(double mass)
{
    double dm = mass - _mass_empty;
    for(auto i = 0u; i < 3; i++)
    {
        _MOI[i] = _MOI_empty[i] + _MOI_delta[i]*dm;
        _COM[i] = _COM_empty[i] + _COM_delta[i]*dm;
    }
    _MOI[3] = _MOI_empty[3] + _MOI_delta[3]*dm;
}

Eigen::Vector3d VehicleBase::acceleration_in_ecef() const 
{
    Eigen::Vector3d ecef_acceleration = _body_frame_ecef.transpose()*(_body_force/_state.body.mass);
    ecef_acceleration += Earth::get_fictional_forces(_state.body.position, _state.body.velocity) + Earth::get_J2_gravity(_state.body.position);
    return ecef_acceleration;
}

Eigen::Vector3d VehicleBase::acceleration_in_body() const 
{
    Eigen::Vector3d angular_accleration;
    functions::angular_acceleration_from_torque_plane_symmetry(_body_moment.data(), 
        _state.body.angular_velocity.data(), _inertia.get_MOI().data(), angular_accleration.data());
    return angular_accleration;
}

void set_acceleration(const Eigen::Vector3d& body_force, const Eigen::Matrix3d& orientation, double mass, 
    const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, double* acc)
{
    const Eigen::Vector3d accel = orientation.transpose()*(body_force/mass);
    const Eigen::Vector3d gravity = Earth::get_J2_gravity(position);
    const double x_fic = Earth::EARTH_ROTATION_RATE*(Earth::EARTH_ROTATION_RATE*position[0] - 2.0*velocity[1]);
    const double y_fic = Earth::EARTH_ROTATION_RATE*(Earth::EARTH_ROTATION_RATE*position[1] + 2.0*velocity[0]);
    acc[0] = accel[0] + gravity[0] - x_fic;
    acc[1] = accel[1] + gravity[1] - y_fic;
    acc[2] = accel[2] + gravity[2];
}

void VehicleBase::set_dx(std::array<double, 14>& dx)
{
    Eigen::Vector3d x = _body_frame_ecef.row(0);
    memcpy(&dx[0], _state.body.velocity.data(), 3*sizeof(double));
    set_acceleration(_body_force, _body_frame_ecef, _state.body.mass, _state.body.position, _state.body.velocity, &dx[3]);
    functions::quaternion_orientation_rate(&_state.x[6],_state.body.angular_velocity.data(), &dx[6]);
    functions::angular_acceleration_from_torque_plane_symmetry(_body_moment.data(), 
        _state.body.angular_velocity.data(), _inertia.get_MOI().data(), &dx[10]);
    dx[13] = _mass_rate;
}

void VehicleBase::update(double time)
{
    _inertia.set_mass(_state.body.mass);
    this->update_environment();
    this->update_control(time);
    this->update_body_forces(time);
}

void VehicleBase::update_environment()
{
    // Get frames, rotation matrices, and other reference quantities
    // Body Frame
    _body_frame_ecef = _state.body.orientation.conjugate().toRotationMatrix();
    // Geodetic Frames
    Earth::getGeodeticFrame(_state.body.position, _lla, _ENU2ECEF.data());
    // Get Air 
    _atmosphere.set_air(_lla.altitude, _air);
    // Get the aero quantities used by aerodynamics
    _aero.update(_air, _state.body.velocity, _body_frame_ecef);
}



VehicleBase::operator bool() const
{
    return _lla.altitude > -1.0 & _lla.altitude < 1e8;
}

void Vehicle::update_control(double time)
{
    gnc->update(*this, time);
}

void Vehicle::update_body_forces(double time)
{
    aerodynamics->update_forces(_aero);
    propulsion.update_forces(_air, _aero, _inertia.get_COM(), time);

    _body_force = aerodynamics->get_body_force() + propulsion.get_body_force();
    _body_moment = aerodynamics->get_body_moment() + propulsion.get_body_moment();

    _mass_rate = -propulsion.thruster->get_mass_rate();
}

void BasicVehicle::update_control(double)
{
    double weight = _state.body.mass*GForce::G;
    double lift = weight + _thruster.get_thrust() - _aerodynamics.get_drag();
    _aerodynamics.set_lift(lift);
}

void BasicVehicle::update_body_forces(double)
{
    _aerodynamics.update_forces(_aero); 

    _body_force = _aerodynamics.get_body_force();
    _body_force[0] += _thruster.get_thrust();

    _mass_rate = -_thruster.get_mass_rate();
}

RamjetVehicle::RamjetVehicle(const InertialProperties& I, const Atmosphere& atmosphere, 
    std::unique_ptr<Ramjet> ramjet, 
    const AerodynamicBasicCoefficients::Coef& coef,
    const AltitudeControl& control) : 
        VehicleBase(I,atmosphere),
        _ramjet(std::move(ramjet)),
        _aerodynamics(coef),
        _control(control)
{}

double AltitudeControl::update_elevator(double time, double altitude, double altitude_rate,
    double airspeed, double acceleration, double mass, double AoA, const double pitch_rate)
{
    const auto altitude_rate_err = std::max(_min_altitude_rate - altitude_rate, 0.0);
    double desired_pitch_rate = altitude_rate_err*_K1 + (_acceleration_desired - acceleration)*_K2;

    desired_pitch_rate = std::clamp(desired_pitch_rate, _max_pitch_down, _max_pitch_up);

    const auto pitch_rate_err = desired_pitch_rate - pitch_rate;
    
    const auto elevator_pitching = pitch_rate_err*_K4;

    const auto alpha_above_max = _max_alpha - AoA;
    const auto alpha_below_min = _min_alpha - AoA; // min is always 0, never have negative lift
    const auto elevator_alpha = (std::min(alpha_above_max, 0.0) + std::max(alpha_below_min, 0.0))*_alpha_k;

    const auto elevator_deflection = elevator_pitching + elevator_alpha;
    
    #ifdef DEBUG
        double pitch = asin(2*(_state.body.orientation.x()*_state.body.orientation.z() 
                 + _state.body.orientation.y()*_state.body.orientation.w()))*57.295779513082;
        std::cout << "current pitch: " << pitch << " ";
        std::cout << "current climb right: " << current_climb_rate << " ";
        std::cout << "desired: " << desired_climb_rate << " ";
        std::cout << "elevator: " << elevator_deflection << " ";
    #endif

    return elevator_deflection;
}

void RamjetVehicle::update_control(double time)
{
    const auto dt = time - _old_time;
    if(fabs(dt) < 1e-3) 
    {
        return;
    }
    const auto inv_dt = 1.0/dt;
    const auto altitude = _lla.altitude;
    const auto airspeed = _aero.airspeed;
    const auto altitude_rate = (altitude - _old_altitude)*inv_dt;
    const auto airspeed_rate = (airspeed - _old_airspeed)*inv_dt;
    const auto mass = _state.body.mass;
    const auto AoA = _aero.alpha_angle;
    const auto pitch_rate = _state.body.angular_velocity.y();

    _old_time = time;
    _old_airspeed = airspeed;
    _old_altitude = altitude;

    const auto elevator = _control.update_elevator(time, altitude, altitude_rate, airspeed, 
        airspeed_rate, mass, AoA, pitch_rate);
    _aerodynamics.set_elevator(elevator);
}

void RamjetVehicle::update_body_forces(double time)
{
    #ifdef DEBUG
        std::cout << "\n*** Time: " << time << "***\n";
    #endif

    _aerodynamics.update_forces(_aero); 

    _body_force = _aerodynamics.get_body_force();
    _body_moment = _aerodynamics.get_body_moment();

    _ramjet->update_thrust(_air, _aero);

    const double thrust = _ramjet->get_thrust();
    if(thrust > _state.body.mass) 
    {
        _body_force[0] += thrust;
        _mass_rate = -_ramjet->get_mass_rate();
    } 
    else 
    {
        _mass_rate = 0;
    }
    

    #ifdef DEBUG
        std::cout << "Thrust: " << _ramjet.get_thrust() << " ";
        std::cout << "ISP: " << _ramjet.get_thrust()/(_ramjet.get_mass_rate()*9.806) << " ";
        std::cout << "Dyn Pres: " <<_aero.dynamic_pressure;
        std::cout << std::endl;
    #endif
}

std::unique_ptr<double[]> RamjetVehicle::get_data() const 
{
    std::unique_ptr<double[]> output = std::make_unique<double[]>(NUM_DATA);
    output[0] = _aerodynamics.get_elevator();
    output[1] = _ramjet->get_thrust();
    return output;
}
