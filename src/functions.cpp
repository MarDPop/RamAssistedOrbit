#include "functions.hpp"

#include <cmath>

namespace functions
{

void quaternion_orientation_rate(const double* q, 
    const double* angular_velocity_body, double* q_dot)
{
    constexpr unsigned int X = 0;
    constexpr unsigned int Y = 1;
    constexpr unsigned int Z = 2;
    constexpr unsigned int W = 3;
    q_dot[X] = (angular_velocity_body[Z]*q[Y] - angular_velocity_body[Y]*q[Z] + angular_velocity_body[X]*q[W])*0.5;
    q_dot[Y] = (-angular_velocity_body[Z]*q[X] + angular_velocity_body[X]*q[Z] + angular_velocity_body[Y]*q[W])*0.5;
    q_dot[Z] = (angular_velocity_body[Y]*q[X] - angular_velocity_body[X]*q[Y] + angular_velocity_body[Z]*q[W])*0.5;
    q_dot[W] = (-angular_velocity_body[X]*q[X] - angular_velocity_body[Y]*q[Y] - angular_velocity_body[Z]*q[Z])*0.5; // scalar
}

void quaternion_orientation_rate(const Eigen::Quaterniond& q, 
    const Eigen::Vector3d& angular_velocity_body, Eigen::Quaterniond& q_dot)
{
    Eigen::Quaterniond w(0.0, 0.5*angular_velocity_body[0], 0.5*angular_velocity_body[1],
        0.5*angular_velocity_body[2]);
    q_dot = q*w;
}

void angular_acceleration_from_torque_principal_axis(const double* torque, const double* angular_velocity,
        const double* inertia, double* angular_acceleration)
{
    angular_acceleration[0] = (torque[0] - (inertia[2] - inertia[1])*angular_velocity[2]*angular_velocity[1])/inertia[0];
    angular_acceleration[1] = (torque[1] - (inertia[0] - inertia[2])*angular_velocity[2]*angular_velocity[0])/inertia[1];
    angular_acceleration[2] = (torque[2] - (inertia[1] - inertia[0])*angular_velocity[1]*angular_velocity[0])/inertia[2];
}

void angular_acceleration_from_torque_plane_symmetry(const double* torque, const double* angular_velocity,
        const double* inertia, double* angular_acceleration)
{
    double xTerm = angular_velocity[1]*angular_velocity[0]*inertia[3];
    double yTerm = (angular_velocity[2]*angular_velocity[2] - angular_velocity[0]*angular_velocity[0])*inertia[3];
    double zTerm = -angular_velocity[1]*angular_velocity[2]*inertia[3];

    double L = torque[0] + (inertia[1] - inertia[2])*angular_velocity[2]*angular_velocity[1] + xTerm;
    double M = torque[1] + (inertia[2] - inertia[0])*angular_velocity[2]*angular_velocity[0] + yTerm;
    double N = torque[2] + (inertia[0] - inertia[1])*angular_velocity[1]*angular_velocity[0] + zTerm;

    angular_acceleration[1] = M/inertia[1];

    double invDet = 1.0/(inertia[0]*inertia[2] - inertia[3]*inertia[3]);

    angular_acceleration[0] = (L*inertia[2] + N*inertia[3])*invDet;
    angular_acceleration[2] = (L*inertia[3] + N*inertia[0])*invDet;
}

// angle Y +/- 90 deg
void ZY_rotation(double yaw, double pitch, double* rotm)
{
    // for intrinsic rotation order is Rz*Ry
    // Eigen Matrix3d stored column major by default
    rotm[2] = -sin(pitch); // xz
    rotm[3] = -sin(yaw); // yx
    rotm[8] = sqrt(1.0 - rotm[2]*rotm[2]); // zz 
    rotm[4] = sqrt(1.0 - rotm[3]*rotm[3]); // yy
    rotm[5] = 0.0; // yz

    rotm[0] = rotm[4]*rotm[8]; // xx
    rotm[1] = -rotm[3]*rotm[8]; // xy
    rotm[6] = -rotm[4]*rotm[2]; // zx
    rotm[7] = rotm[3]*rotm[2]; // zy
    
}

} // end namespace