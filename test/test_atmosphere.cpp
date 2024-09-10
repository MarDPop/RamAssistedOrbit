#include "doctest.h"

#include "atmosphere.hpp"

TEST_CASE("Area Ratio")
{
    double mach_test = 1.2;
    double area_ratio = Air::isentropic_area_ratio(mach_test);
    REQUIRE(area_ratio == doctest::Approx(0.9704594567513530).epsilon(1e-6));
    double mach = Air::supersonic_mach_area_ratio(area_ratio, 1.4);
    REQUIRE(mach  == doctest::Approx(mach_test).epsilon(1e-4));

    mach_test = 4.0;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    REQUIRE(area_ratio == doctest::Approx(0.0932944606413994169).epsilon(1e-6));
    mach = Air::supersonic_mach_area_ratio(area_ratio, 1.4);
    REQUIRE(mach  == doctest::Approx(mach_test).epsilon(1e-4));

    mach_test = 7.5;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    REQUIRE(area_ratio == doctest::Approx(0.0070501236933626814).epsilon(1e-6));
    mach = Air::supersonic_mach_area_ratio(area_ratio, 1.4);
    REQUIRE(mach == doctest::Approx(mach_test).epsilon(1e-4));

    mach_test = 0.8;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    REQUIRE(area_ratio == doctest::Approx(0.96317771592).epsilon(1e-6));
    mach = Air::subsonic_mach_area_ratio(area_ratio, 1.4);
    REQUIRE(mach == doctest::Approx(mach_test).epsilon(1e-4));

    mach_test = 0.2;
    area_ratio = Air::isentropic_area_ratio(mach_test);
    REQUIRE(area_ratio == doctest::Approx(0.33743656192).epsilon(1e-6));
    mach = Air::subsonic_mach_area_ratio(area_ratio, 1.4);
    REQUIRE(mach == doctest::Approx(mach_test).epsilon(1e-4));
}