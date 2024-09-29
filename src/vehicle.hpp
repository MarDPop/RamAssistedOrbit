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

    const Atmosphere& _atmosphere;

public:

    VehicleBase(const InertialProperties& I, const Atmosphere& atmosphere) : 
        _inertia(I), _atmosphere(atmosphere) {}

    virtual ~VehicleBase() {}

    [[nodiscard]] Eigen::Vector3d acceleration_in_ecef() const;

    [[nodiscard]] Eigen::Vector3d acceleration_in_body() const;

    [[nodiscard]] const Air& get_air() const
    {
        return _air;
    }

    [[nodiscard]] const AeroQuantities& get_aero() const
    {
        return _aero;
    }

    [[nodiscard]] const Earth::Geodetic& get_lla() const
    {
        return _lla;
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

class AltitudeControl 
{
    const double _K1;

    const double _K2;

    const double _K3;

    const double _K4;

    const double _K5;

    const double _max_alpha;

    const double _min_alpha;

    const double _alpha_k;

    const double _cruise_altitude;

    const double _max_pitch_up;

    const double _max_pitch_down;

    const double _cruise_climb_rate;

    const double _climb_min_rate;

    bool _cruising = false;

public:

    AltitudeControl(double K1, double K2, double K3, double K4, double K5,
        double max_alpha, double min_alpha, double alpha_k,
        double cruise_altitude, double max_pitch_up = 0.1, double max_pitch_down = -0.1,
        double cruise_climb_rate = 1.0, double min_altitude_rate = 20.0) : 
            _K1(K1), _K2(K2), _K3(K3), _K4(K4), _K5(K5), _max_alpha(max_alpha), _min_alpha(min_alpha), _alpha_k(alpha_k),
            _cruise_altitude(cruise_altitude), _max_pitch_up(max_pitch_up),
            _max_pitch_down(max_pitch_down), _cruise_climb_rate(cruise_climb_rate), _climb_min_rate(min_altitude_rate) {}

    AltitudeControl(const AltitudeControl& other) = default;
    AltitudeControl(AltitudeControl&& other) noexcept = default;

    void reset() noexcept 
    {
        _cruising = false;
    }

    double update_elevator(double time, double altitude, double altitude_rate,
        double airspeed, double acceleration, double AoA, double pitch, double pitch_rate, double dynamic_pressure);

};

class RamjetVehicle final : public virtual VehicleBase
{
    std::unique_ptr<Ramjet> _ramjet;

    AerodynamicBasicCoefficients _aerodynamics;

    AltitudeControl _control;

    double _old_time = 0.0;

    double _old_altitude = 0.0;

    double _old_airspeed = 0.0;

public:

    static constexpr unsigned NUM_DATA = 2u;

    RamjetVehicle(const InertialProperties& I, const Atmosphere& atmosphere, 
        std::unique_ptr<Ramjet> ramjet, 
        const AerodynamicBasicCoefficients::Coef& coef,
        const AltitudeControl& control);

    void initNav(double time, double altitude, double airspeed)
    {
        _old_time = time;
        _old_altitude = altitude;
        _old_airspeed = airspeed;
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