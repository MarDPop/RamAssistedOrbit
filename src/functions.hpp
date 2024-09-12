#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Eigen/Dense>

namespace functions
{
    

    /**
     * 
     * Note: all expressed in body
     */
    void angular_acceleration_from_torque_principal_axis(const double* torque, const double* angular_velocity,
        const double* inertia, double* angular_acceleration);

    Eigen::Vector3d angular_acceleration_from_torque_principal_axis(const Eigen::Vector3d& torque, const Eigen::Vector3d& angular_velocity,
        const std::array<double, 3>& inertia);

    void angular_acceleration_from_torque_plane_symmetry(const double* torque, const double* angular_velocity,
        const double* inertia, double* angular_acceleration);

    Eigen::Vector3d angular_acceleration_from_torque_plane_symmetry(const Eigen::Vector3d& torque, const Eigen::Vector3d& angular_velocity,
        const std::array<double, 4>& inertia);

    void normalize_quat(double* quaternion);

    void quaternion_orientation_rate(const double* current_inertial_orientation_quaternion, 
        const double* angular_velocity_body, double* inertial_quaternion_rate);

    void quaternion_orientation_rate(const Eigen::Quaterniond& q, 
        const Eigen::Vector3d& angular_velocity_body, Eigen::Quaterniond& q_dot);

    /**
     * 
     * Note: yaw followed by pitch intrinsic, rotm stored in col major
     */
    void ZY_rotation(double yaw, double pitch, double* rotm);
}

#endif