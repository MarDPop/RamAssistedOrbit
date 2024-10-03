#ifndef VEHICLE_H
#define VEHICLE_H

#include "aerodynamics.hpp"
#include "earth.hpp"
#include "functions.hpp"
#include "gnc.hpp"
#include "ode.hpp"
#include "propulsion.hpp"
#include "state_6dof.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Dense"

#include <array>
#include <exception>
#include <memory>
#include <string>

class InertialProperties
{
    Eigen::Vector3d _COM;

    std::array<double, 4> _MOI;

    std::array<double, 4> _MOI_delta;

    Eigen::Vector3d _COM_delta;

public:

    const double mass_empty;

    const std::array<double, 4> MOI_empty;

    const Eigen::Vector3d COM_empty;

    static constexpr double IXX = 0;
    static constexpr double IYY = 1;
    static constexpr double IZZ = 2;

    InertialProperties(double mass, std::array<double, 4> MOI, Eigen::Vector3d COM);

    InertialProperties(double mass_empty, std::array<double, 4> MOI_empty, Eigen::Vector3d COM_empty,
        double mass_full, const std::array<double, 4>& MOI_full, const Eigen::Vector3d& COM_full);

    [[nodiscard]] const Eigen::Vector3d& get_COM() const
    {
        return _COM;
    }

    [[nodiscard]] const std::array<double, 4>& get_MOI() const
    {
        return _MOI;
    }

    void set_mass(double mass);
};

class VehicleBase : public virtual Dynamics<14>
{
protected: 

    friend class Navigation;

    friend class Vehicle_ODE;

    State_6DOF _state{};

    Eigen::Matrix3d _body_frame_ecef = Eigen::Matrix3d::Identity();

    InertialProperties _inertia;

    Eigen::Vector3d _body_force = {0.0,0.0,0.0};

    Eigen::Vector3d _body_moment = {0.0,0.0,0.0};

    double _mass_rate = 0.0;

    Earth::Geodetic _lla = {0.0, 0.0, 0.0};

    Eigen::Matrix3d _ENU2ECEF = Eigen::Matrix3d::Identity();

    Air _air{};

    AeroQuantities _aero{};

    const Atmosphere* _atmosphere;

public:

    VehicleBase(const InertialProperties& I, const Atmosphere& atmosphere) : 
        _inertia(I), _atmosphere(&atmosphere) {}

    virtual ~VehicleBase() {}

    [[nodiscard]] Eigen::Vector3d acceleration_in_ecef() const;

    [[nodiscard]] Eigen::Vector3d acceleration_in_body() const;

    [[nodiscard]] const AeroQuantities& get_aero() const
    {
        return _aero;
    }

    [[nodiscard]] const Air& get_air() const
    {
        return _air;
    }

    [[nodiscard]] const Eigen::Vector3d& get_body_force() const
    {
        return _body_force;
    }

    [[nodiscard]] const Eigen::Vector3d& get_body_moment() const
    {
        return _body_moment;
    }

    [[nodiscard]] const Eigen::Matrix3d& get_body_frame_ecef() const
    {
        return _body_frame_ecef;
    }

    [[nodiscard]] const Eigen::Matrix3d& get_ENU2ECEF() const
    {
        return _ENU2ECEF;
    }
    
    [[nodiscard]] const InertialProperties& get_inertia() const
    {
        return _inertia;
    }
    
    [[nodiscard]] const Earth::Geodetic& get_lla() const
    {
        return _lla;
    }

    [[nodiscard]] double get_mass_rate() const
    {
        return _mass_rate;
    }

    [[nodiscard]] const State_6DOF& get_state() const
    {
        return _state;
    }

    [[nodiscard]] std::vector<double> get_other_values() const 
    {
        return std::vector<double>{};
    }

    void update_environment();

    void set_atmosphere(const Atmosphere* atmosphere)
    {
        if(atmosphere == nullptr)
        {
            throw std::invalid_argument("Atmosphere cannot be null.");
            return;
        }
        _atmosphere = atmosphere;
    }

    void set_dx(std::array<double, 14>& dx);

    void update(double time);

    virtual void update_control([[maybe_unused]] double time) {}

    virtual void update_body_forces([[maybe_unused]] double time) {}

    void operator()(std::array<double, 14>& x, const double t, 
        std::array<double, 14>& dx) override final
    {
        // Set state
        functions::normalize_quat(&x[6]);
        _state.x = x;
        this->update(t);
        this->set_dx(dx);
    }

    operator bool() const override;

};

class Vehicle final : public virtual VehicleBase
{

public:

    std::unique_ptr<Aerodynamics> const aerodynamics;

    std::unique_ptr<VehicleGNC> const gnc;

    PropulsionSystem propulsion;

    Vehicle(const InertialProperties& I, const Atmosphere& atmosphere, 
        std::unique_ptr<Aerodynamics> aerodynamics_, std::unique_ptr<VehicleGNC> gnc_, std::unique_ptr<Thruster> thruster_, 
        std::unique_ptr<ThrusterControl> control_) : VehicleBase(I,atmosphere), aerodynamics(std::move(aerodynamics_)), 
            gnc(std::move(gnc_)), propulsion(std::move(thruster_), std::move(control_))
    {
        if(!aerodynamics || !gnc)
        {
            throw std::runtime_error("Invalid aerodynamics");
        }
    }

    void update_control(double time) override;

    void update_body_forces(double time) override;

};

class BasicVehicle final : public virtual VehicleBase 
{
    Thruster _thruster;

    Lift2DragAerodynamics _aerodynamics;

public:

    BasicVehicle(const InertialProperties& I, const Atmosphere& atmosphere, 
        const Thruster& thruster, const double L2D) 
            : VehicleBase(I, atmosphere), _thruster(thruster), _aerodynamics(L2D) {}

    void update_control(double time) override;

    void update_body_forces(double time) override;
};

class PitchGuidance
{
protected:

    double _desired_pitch;

public:

    double get_desired_pitch() const
    {
        return _desired_pitch;
    }
    static double opt_alpha(double CL_alpha, double K, double CD_0);

    void set_desired_pitch(double pitch)
    {
        _desired_pitch = pitch;
    }
    
    virtual void update_climb_navigation(const VehicleBase& vehicle, double time){}

};

class PitchControl
{
    double _pitch_commanded;

public:

    const double K;

    const double D;

    PitchControl(double K_, double D_) : K(K_), D(D_) {}

    void set_pitch_commanded(double pitch) 
    {
        _pitch_commanded = pitch;
    }

    double get_elevator(double pitch, double pitch_rate, double dynamic_pressure) const;
};

class AltitudeRateGuidance : public virtual PitchGuidance
{
    const double _accel_K;

    const double _alpha_accel_K;

    const double _altitude_K;

    const double _max_alpha;

    const double _min_alpha;

    const double _opt_alpha;

    const double _alpha_k;

    const double _cruise_altitude;

    const double _cruise_climb_rate;

    const double _climb_min_rate;

    const double _max_pitch_offset;

    double _pitch = 0.0;

    double _pitch_rate = 0.0;

    double _old_time = 0.0;

    double _old_altitude = 0.0;

    double _old_airspeed = 0.0;

    bool _cruising = false;

public:

    AltitudeRateGuidance(double K1, double K2, double K3,
        double max_alpha, double min_alpha, double opt_alpha, double alpha_k,
        double cruise_altitude, double max_pitch_offset,
             double cruise_climb_rate = 1.0, double min_climb_rate = 20.0) : 
            _accel_K(K1), _alpha_accel_K(K2), _altitude_K(K3), _max_alpha(max_alpha), _min_alpha(min_alpha), _opt_alpha(opt_alpha),
            _alpha_k(alpha_k), _cruise_altitude(cruise_altitude), _max_pitch_offset(max_pitch_offset), 
            _cruise_climb_rate(cruise_climb_rate), _climb_min_rate(min_climb_rate) {}

    AltitudeRateGuidance(const AltitudeRateGuidance& other) = default;
    AltitudeRateGuidance(AltitudeRateGuidance&& other) noexcept = default;

    double get_pitch() const
    {
        return _pitch;
    }

    double get_pitch_rate() const
    {
        return _pitch_rate;
    }

    void init(const VehicleBase& vehicle, double time);

    void reset() noexcept 
    {
        _cruising = false;
    }

    void update_climb_navigation(const VehicleBase& vehicle, double time) override;

};

class RamjetVehicle final : public virtual VehicleBase
{
    std::unique_ptr<Ramjet> _ramjet;

    AerodynamicBasicCoefficients _aerodynamics;

    AltitudeRateGuidance _guidance;

    PitchControl _control;

public:

    static constexpr unsigned NUM_DATA = 2u;

    RamjetVehicle(const InertialProperties& I, const Atmosphere& atmosphere, 
        std::unique_ptr<Ramjet> ramjet, 
        const AerodynamicBasicCoefficients::Coef& coef,
        const AltitudeRateGuidance& guidance,
        const PitchControl& control);

    void initNav(double time = 0.0)
    {
        _guidance.init(*this, time);
    }

    void update_control(double time) override;

    void update_body_forces(double time) override;

    unsigned num_data() const override 
    {
        return NUM_DATA;
    }

    std::unique_ptr<double[]> get_data() const override;

};

#endif