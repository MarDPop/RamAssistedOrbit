#include "test.hpp"

#include "earth.hpp"
#include "functions.hpp"
#include "propulsion.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Dense"

#include <iostream>
#include <memory>
#include <vector>

#define DEG2RAD 1.7453292519943295769e-2

bool TestMath::testZYRotation()
{
    double zAngle = 0.2;
    double yAngle = 0.1;
    Eigen::Matrix3d Zrotation = Eigen::AngleAxisd(zAngle, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    Eigen::Matrix3d Yrotation = Eigen::AngleAxisd(yAngle, Eigen::Vector3d::UnitY()).toRotationMatrix();

    Eigen::Matrix3d ZYrotation = Zrotation*Yrotation;

    Eigen::Matrix3d myZY;
    functions::ZY_rotation(zAngle, yAngle, myZY.data());

    Eigen::Matrix3d diff = ZYrotation - myZY;

    for(auto i = 0u; i < 9; i++)
    {
        assertZero(diff.data()[i], 1e-6);
    }
    return true;
}

bool TestEarth::testLLA()
{
    std::array<double, 3> expectedLLA = {30.0*DEG2RAD, 20.0*DEG2RAD, 550.0};

    Eigen::Vector3d ecef(5192994.214, 2998176.607, 2167884.90);
    const auto lla = Earth::ECEF2LLA(ecef);

    assertClose(expectedLLA[0], lla.longitude, 1e-7);
    assertClose(expectedLLA[1], lla.latitude, 1e-7);
    assertClose(expectedLLA[2], lla.altitude, 1e-2);

    Eigen::Vector3d expectedUP(cos(expectedLLA[0])*cos(expectedLLA[1]),sin(expectedLLA[0])*cos(expectedLLA[1]), sin(expectedLLA[1]));

    auto up = Earth::getUp(ecef);

    assertArrayClose(expectedUP.data(), up.data(), 3, 1e-6);

    return true;
}

bool TestPhysics::test_area_ratio()
{
    bool flag = true;
    double mach_test = 1.2;
    double area_ratio = Air::isentropic_area_ratio(mach_test);
    flag &= assertClose(area_ratio, 0.9704594567513530, 1e-6);
    double mach = Air::supersonic_mach_area_ratio(area_ratio, 1.4);
    flag &= assertClose(mach, mach_test, 1e-4);

    mach_test = 4.0;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    flag &= assertClose(area_ratio, 0.0932944606413994169, 1e-6);
    mach = Air::supersonic_mach_area_ratio(area_ratio, 1.4);
    flag &= assertClose(mach, mach_test, 1e-4);

    mach_test = 7.5;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    flag &= assertClose(area_ratio, 0.0070501236933626814, 1e-6);
    mach = Air::supersonic_mach_area_ratio(area_ratio, 1.4);
    flag &= assertClose(mach, mach_test, 1e-4);

    mach_test = 0.8;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    flag &= assertClose(area_ratio, 0.96317771592, 1e-6);
    mach = Air::subsonic_mach_area_ratio(area_ratio, 1.4);
    flag &= assertClose(mach, mach_test, 1e-4);

    mach_test = 0.2;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    flag &= assertClose(area_ratio, 0.33743656192, 1e-6);
    mach = Air::subsonic_mach_area_ratio(area_ratio, 1.4);
    flag &= assertClose(mach, mach_test, 1e-4);

    return flag;
}

bool TestPhysics::test_aero_quantities()
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

    flag &= assertClose(aero.airspeed, airspeed1, 1e-6);
    flag &= assertClose(aero.mach, 0.5, 1e-6);
    flag &= assertClose(aero.dynamic_pressure, 5000, 1e-6);
    flag &= assertClose(aero.alpha_angle, AoA1, 1e-6);
    flag &= assertClose(aero.beta_angle, 0, 1e-6);

    return flag;
}

bool TestRamjet::testRamjet()
{
    double cruise_mach = 4.5;
    double L2D = 4*(cruise_mach + 3)/cruise_mach;
    auto ramjet = RamjetVariableInlet::create(0.6, 20000, 4.0, 4.5, L2D, 10000, 0.25);

    std::cout << ramjet.get_mass_rate() << " " << ramjet.get_thrust() << std::endl;

    return true;
}

void Tests::runAll()
{
    std::vector<std::unique_ptr<Test>> tests;
    TestMath test1;
    test1.run();

    TestEarth test2;
    test2.run();

    TestPhysics test3;
    test3.run();

    TestRamjet test4;
    test4.run();
}