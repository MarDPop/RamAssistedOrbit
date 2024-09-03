#ifndef VEHICLE_H
#define VEHICLE_H

#include "aerodynamics.hpp"
#include "earth.hpp"
#include "functions.hpp"
#include "gnc.hpp"
#include "ode.hpp"
#include "propulsion.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Dense"

#include <array>
#include <exception>
#include <memory>
#include <string>

union State
{
    std::array<double, 14> x;
    struct
    {
        Eigen::Vector3d position; // in ECEF meters
        Eigen::Vector3d velocity; // in ECEF meters / sec
        Eigen::Quaterniond orientation; // remember stored x,y,z,w 
        // represents rotation matrix [x axis forward, y to right, z down] (row major, ie rotation from ecef to body frame)
        Eigen::Vector3d angular_velocity; // remember this is expressed in body frame in rad / sec
        double mass; // in kg
    } body; 

    State() {}
    State(const State& state) : x(state.x) {}

    [[nodiscard]] std::string to_json_string() const;
};

class InertialProperties
{
    Eigen::Vector3d _COM;

    std::array<double, 4> _MOI;

    const double _mass_empty;

    const std::array<double, 4> _MOI_empty;

    const Eigen::Vector3d _COM_empty;

    std::array<double, 4> _MOI_delta;

    Eigen::Vector3d _COM_delta;

public:

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

    State _state;

    Eigen::Matrix3d _body_frame_ecef;

    InertialProperties _inertia;

    Eigen::Vector3d _body_force;

    Eigen::Vector3d _body_moment;

    double _mass_rate = 0.0;

    Earth::Geodetic _lla;

    Eigen::Matrix3d _ENU2ECEF = Eigen::Matrix3d::Identity();

    Air _air;

    AeroQuantities _aero;

    const Atmosphere& _atmosphere;

public:

    VehicleBase(const InertialProperties& I, const Atmosphere& atmosphere) : 
        _inertia(I), _body_force(0,0,0), _body_moment(0,0,0), _atmosphere(atmosphere) {}

    virtual ~VehicleBase() {}

    [[nodiscard]] const Air& get_air() const
    {
        return _air;
    }

    [[nodiscard]] const AeroQuantities& get_aero() const
    {
        return _aero;
    }

    [[nodiscard]] double get_altitude() const
    {
        return _lla.altitude;
    }

    void update_environment();

    void set_dx(std::array<double, 14>& dx);

    virtual void update_control([[maybe_unused]] double time) {}

    virtual void update_body_forces([[maybe_unused]] double time) {}

    void operator()(const std::array<double, 14>& x, const double t, 
        std::array<double, 14>& dx) override final;

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

class RamjetVehicle final : public virtual VehicleBase
{
    RamjetFixedInlet _ramjet;

    AerodynamicBasicCoefficients _aerodynamics;

    double _K1;

    double _C1;

    double _max_alpha;

    double _min_alpha;

    double _alpha_k;

    double _cruise_altitude;

    double _cruise_mach;

public:

    static constexpr unsigned NUM_DATA = 2u;

    RamjetVehicle(const InertialProperties& I, const Atmosphere& atmosphere, 
        const RamjetFixedInlet& ramjet, const AerodynamicBasicCoefficients::Coef& coef);
    
    void set_control_values(double K1, double C1, double max_alpha, double min_alpha, double alpha_k,
        double cruise_altitude, double cruise_mach)
    {
        _K1 = K1;
        _C1 = C1;
        _max_alpha = max_alpha;
        _min_alpha = min_alpha;
        _alpha_k = alpha_k;
        _cruise_altitude = cruise_altitude;
        _cruise_mach = cruise_mach;
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