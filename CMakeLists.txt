cmake_minimum_required(VERSION 3.12)

project(AirBreathingAssist)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/lib/Eigen)

file(GLOB SRC_FILES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}"/src/*.hpp" ${PROJECT_SOURCE_DIR}"/src/*.cpp")

add_executable( sim.exe ${SRC_FILES})