#ifndef FUNCTIONS_H
#define FUNCTIONS_H

namespace functions
{
    /**
     * 
     */
    void quaternion_orientation_rate(const double* current_inertial_orientation_quaternion, 
        const double* angular_velocity_body, double* inertial_quaternion_rate);

    /**
     * 
     * Note: all expressed in body
     */
    void angular_acceleration_from_torque_principal_axis(const double* torque, const double* angular_velocity,
        const double* inertia, double* angular_acceleration);

    /**
     * 
     * Note: all expressed in body
     */
    void angular_acceleration_from_torque_plane_symmetry(const double* torque, const double* angular_velocity,
        const double* inertia, double* angular_acceleration);

    /**
     * 
     * Note: yaw followed by pitch intrinsic, rotm stored in col major
     */
    void ZY_rotation(double yaw, double pitch, double* rotm);
}

#endif