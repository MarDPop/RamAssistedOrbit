#include "aerodynamics.hpp"

#include <cmath>

void AeroQuantities::update(const Air& air, const Eigen::Vector3d& velocity_ecef, const Eigen::Matrix3d& orientation)
{
    // assert(air_body_vector.x() > 0.0);

    // remember x axis forward, z axis down

    // no wind
    airspeed = velocity_ecef.norm();
    mach = airspeed*air.inv_sound_speed;
    dynamic_pressure = 0.5*air.density*airspeed*airspeed; // or air.gamma*0.5*air.pressure*mach*mach

    // relative air vector (multiplied by -airspeed to do) (no wind)
    air_body_vector = (orientation*velocity_ecef)/-airspeed;

    // relative air angles
    alpha_angle = atan(air_body_vector.z()/air_body_vector.x());
    beta_angle = atan(air_body_vector.y()/air_body_vector.x());

    // Lift vector is v cross y
    auto lift_normalization = 1.0/sqrt(air_body_vector.x()*air_body_vector.x() + air_body_vector.z()*air_body_vector.z());
    body_lift_vector(0) = -air_body_vector.z()*lift_normalization;
    body_lift_vector(2) = air_body_vector.x()*lift_normalization;

    reynolds_number_per_meter = air.density*airspeed/air.dynamic_viscosity;
}

void BallisticAerodynamics::update_forces(const AeroQuantities& aero)
{
    _body_force = aero.air_body_vector*(_CDA*aero.dynamic_pressure);
}

void Lift2DragAerodynamics::update_forces(const AeroQuantities& aero)
{
    _body_force = aero.body_lift_vector*_lift - aero.air_body_vector*_drag;
}

void AerodynamicBasicCoefficients::update_forces(const AeroQuantities& aero)
{
    // Base coefficient
    auto CLA = coef_scaled.CL_alpha*(aero.alpha_angle + coef_scaled.alpha0);
    auto CMAl = coef_scaled.CM_alpha*aero.alpha_angle;
    const auto CDA = coef_scaled.CD0 + coef_scaled.K*CLA*CLA;

    if(aero.alpha_angle > coef_scaled.stall_angle)
    {
        double factor = coef_scaled.stall_angle/(aero.alpha_angle*aero.alpha_angle);
        CLA *= factor;
        CMAl *= factor;
    }
    // get elevator contributions
    const auto dCDA = coef_scaled.dCD_dEl*_elevator;
    const auto dCMA = coef_scaled.dCM_dEl*_elevator;

    // get scalar quanties
    const auto drag = (CDA + fabs(dCDA))*aero.dynamic_pressure;
    const auto lift = CLA*aero.dynamic_pressure;
    const auto pitch_moment = (CMAl + dCMA)*aero.dynamic_pressure;

    _body_force = aero.body_lift_vector*lift + aero.air_body_vector*drag;

    // pitch only
    _body_moment(1) = pitch_moment;
}

AerodynamicBasicCoefficients::AerodynamicBasicCoefficients(const Coef& coef) : Aerodynamics(),
    coef_scaled({coef.CD0*coef.A_ref, coef.K/coef.A_ref, coef.alpha0, coef.CL_alpha*coef.A_ref,
    coef.CM_alpha*coef.A_ref*coef.L_ref, coef.dCD_dEl*coef.A_ref, coef.dCM_dEl*coef.A_ref*coef.L_ref,
    coef.stall_angle, coef.max_elevator_deflection, coef.A_ref, coef.L_ref}) {}

const std::array<double, AerodynamicAeroTablePlane::N_MACH> AerodynamicAeroTablePlane::MACH = {0.3, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 3.0, 5.0};

const std::array<double, AerodynamicAeroTablePlane::N_ALPHA> AerodynamicAeroTablePlane::ALPHA = {-0.4, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

const std::array<double, AerodynamicAeroTablePlane::N_BETA> AerodynamicAeroTablePlane::BETA = {0.1, 0.2, 0.4, 0.7};

