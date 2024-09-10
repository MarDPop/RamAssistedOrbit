project(AirBreathingAssist)

set(TEST_FILES test_aerodynamics.cpp test_atmosphere.cpp test_earth.cpp test_functions.cpp)

# The test program
add_executable(tests test_main.cpp ${TEST_FILES})
target_include_directories(tests PUBLIC ${TOP_DIR}/lib/doctest)
target_link_libraries(tests PRIVATE sim_lib)