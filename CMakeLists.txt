cmake_minimum_required(VERSION 3.20)

project(AirBreathingAssist)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

include_directories(${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src)

file(GLOB SRC_FILES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(sim ${SRC_FILES})

# set(CMAKE_CXX_COMPILER g++)
# set(CMAKE_C_COMPILER gcc)
# set(CMAKE_CXX_FLAGS "-static-libgcc -static-libstdc++")
# set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
