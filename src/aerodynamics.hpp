#ifndef AERODYNAMICS_H
#define AERODYNAMICS_H

#include "atmosphere.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Dense>

#include <array>

struct AeroQuantities
{
    double airspeed = 0.0;
    Eigen::Vector3d air_body_vector = Eigen::Vector3d(1,0,0);
    Eigen::Vector3d body_lift_vector = Eigen::Vector3d(0,0,1);
    double mach = 0.0;
    double dynamic_pressure = 0.0;
    double alpha_angle = 0.0; 
    double beta_angle = 0.0;
    double reynolds_number_per_meter = 0.0;

    void update(const Air& air, const Eigen::Vector3d& velocity_ecef, const Eigen::Matrix3d& orientation);
};

class Aerodynamics
{
protected:

    Eigen::Vector3d _body_force;

    Eigen::Vector3d _body_moment;

public:

    Aerodynamics() : _body_force(0,0,0), _body_moment(0,0,0) {}
    virtual ~Aerodynamics() {}

    virtual void update_forces(const AeroQuantities&) {};

    const Eigen::Vector3d& get_body_force() const
    {
        return _body_force;
    }

    const Eigen::Vector3d& get_body_moment() const
    {
        return _body_moment;
    }
};

class BallisticAerodynamics final : public virtual Aerodynamics
{
    const double _CDA;

public:

    BallisticAerodynamics(double CD, double A) : _CDA(-CD*A) {}

    void update_forces(const AeroQuantities& aero) override;

};

class Lift2DragAerodynamics final : public virtual Aerodynamics
{
    double _lift;

    double _drag;

public:

    const double _drag2lift;

    Lift2DragAerodynamics(double L2D) : _drag2lift(1.0/L2D) {}

    void set_lift(double lift)
    {
        _lift = lift;
        _drag = _drag2lift*lift;
    }

    double get_drag() const 
    {
        return _drag;
    }

    double get_lift() const
    {
        return _lift;
    }

    void update_forces(const AeroQuantities& aero) override;

};

class AerodynamicBasicCoefficients final : public virtual Aerodynamics
{

    double _elevator = 0;

public:

    struct Coef
    {
        double CD0;
        double K;
        double alpha0;
        double CL_alpha;
        double CM_alpha;
        double dCD_dEl;
        double dCM_dEl;
        double stall_angle;
        double max_elevator_deflection;
        double A_ref;
        double L_ref;
    };  

    const Coef coef_scaled;

    AerodynamicBasicCoefficients(const Coef& coef);

    double get_elevator() const
    {
        return _elevator;
    }

    void set_elevator(double elevator)
    {
        _elevator = std::clamp(elevator, -coef_scaled.max_elevator_deflection, 
            coef_scaled.max_elevator_deflection);
    }

    void update_forces(const AeroQuantities& aero) override;

};


class AerodynamicAeroTablePlane final : public virtual Aerodynamics
{
    static constexpr unsigned N_MACH = 10;

    static constexpr unsigned N_ALPHA = 12;

    static constexpr unsigned N_BETA = 4;

    std::array<std::array<double, 6>, N_MACH> _CF_nominal;

    std::array<std::array<std::array<double, 6>, N_ALPHA>, N_MACH> _delta_CF_alpha; 

    std::array<std::array<std::array<double, 6>, N_BETA>, N_MACH> _delta_CF_beta;

public:

    static const std::array<double, N_MACH> MACH;

    static const std::array<double, N_ALPHA> ALPHA;

    /**
     * With beta angles, we say it's symmetric and 0 at 0
     */
    static const std::array<double, N_BETA> BETA; 

    void update_forces(const AeroQuantities& aero) override;

};

#endif