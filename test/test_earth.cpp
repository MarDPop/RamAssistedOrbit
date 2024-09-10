#include "doctest.h"

#include "earth.hpp"

#define DEG2RAD 1.7453292519943295769e-2

TEST_CASE("TEST LLA conversion")
{
    std::array<double, 3> expectedLLA = {30.0*DEG2RAD, 20.0*DEG2RAD, 550.0};

    Eigen::Vector3d ecef(5192994.214, 2998176.607, 2167884.90);
    const auto lla = Earth::ECEF2LLA(ecef);

    CHECK(fabs(expectedLLA[0] - lla.longitude) < 1e-7);
    CHECK(fabs(expectedLLA[1] - lla.latitude) < 1e-7);
    CHECK(fabs(expectedLLA[2] - lla.altitude) < 1e-2);

    Eigen::Vector3d expectedUP(cos(expectedLLA[0])*cos(expectedLLA[1]),sin(expectedLLA[0])*cos(expectedLLA[1]), sin(expectedLLA[1]));

    auto up = Earth::getUp(ecef);

    for(int i = 0; i < 3; ++i)
    {
        CHECK(fabs(expectedUP[i] - up[i]) < 1e-6);
    }
}