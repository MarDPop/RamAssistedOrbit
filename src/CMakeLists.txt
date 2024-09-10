project(AirBreathingAssist)

include_directories(${TOP_DIR}/lib)

# All sources that also need to be tested in unit tests go into a static library
set(LIB_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/aerodynamics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/atmosphere.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/earth.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/functions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gnc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/output.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/propulsion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/track_dynamics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vehicle.cpp
)

add_library(sim_lib STATIC ${LIB_SOURCES})
target_include_directories(sim_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(sim_lib INTERFACE ${TOP_DIR}/lib)

# The main program
add_executable(sim main.cpp)
target_link_libraries(sim PRIVATE sim_lib)