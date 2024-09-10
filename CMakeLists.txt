cmake_minimum_required(VERSION 3.20)

project(AirBreathingAssist)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

#SET( CMAKE_CXX_FLAGS_DEBUG  "")
#SET( CMAKE_CXX_FLAGS_RELEASE  "")
# SET( CMAKE_CXX_FLAGS  "-march=native -Wall")

if (MSVC)
    # warning level 4
    # add_compile_options(/W4)
else()
    # additional warnings
    # add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find_package(Eigen3 3.4 REQUIRED)

set(TOP_DIR ${PROJECT_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(test)


# set(CMAKE_CXX_COMPILER g++)
# set(CMAKE_C_COMPILER gcc)
# set(CMAKE_CXX_FLAGS "-static-libgcc -static-libstdc++")
# set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
