#ifndef STATE_6DOF_HPP
#define STATE_6DOF_HPP

#include <array>
#include <string>

#include <Eigen/Dense>

union alignas(16) State_6DOF
{
    std::array<double, 14> x;
    struct alignas(16)
    {
        Eigen::Vector3d position; // in ECEF meters
        Eigen::Vector3d velocity; // in ECEF meters / sec
        Eigen::Quaterniond orientation; // remember stored x,y,z,w 
        // represents rotation matrix [x axis forward, y to right, z down] (row major, ie rotation from ecef to body frame)
        Eigen::Vector3d angular_velocity; // remember this is expressed in body frame in rad / sec
        double mass; // in kg
    } body; 

    State_6DOF() : x{0.0} {}
    State_6DOF(const State_6DOF& state) : x(state.x) {}
    State_6DOF(State_6DOF&& state) : x(std::move(state.x)) {}

    State_6DOF& operator=(const State_6DOF& state) 
    {
        x = state.x;
        return *this;
    }

    State_6DOF& operator=(State_6DOF&& state)
    {
        x = std::move(state.x);
        return *this;
    }
};

#endif