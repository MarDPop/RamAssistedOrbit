#include "doctest.h"

#include "aerodynamics.hpp"

#define DEG2RAD 1.7453292519943295769e-2

TEST_CASE("Test aero quantities calculation with constant air parameters")
{
    bool flag = true;

    Eigen::Matrix3d CS1 = Eigen::Matrix3d::Identity();

    const double AoA1 = 10*DEG2RAD;
    const double airspeed1 = 100;
    Eigen::Vector3d velocity(cos(AoA1)*airspeed1, 0, sin(AoA1)*airspeed1);

    Air air;
    air.density = 1.0;
    air.inv_sound_speed = 1.0/200;
    air.pressure = 100000;

    AeroQuantities aero;

    aero.update(air, velocity, CS1);

    REQUIRE(aero.airspeed == doctest::Approx(airspeed1).epsilon(1e-6));
    REQUIRE(aero.mach  == doctest::Approx( 0.5).epsilon(1e-6));
    REQUIRE(aero.dynamic_pressure == doctest::Approx(5000).epsilon(1e-6));
    REQUIRE(aero.alpha_angle == doctest::Approx(AoA1).epsilon(1e-6));
    REQUIRE(aero.beta_angle == doctest::Approx(0).epsilon(1e-6));
}