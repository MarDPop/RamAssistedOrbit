#ifndef EARTH_H
#define EARTH_H
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <Eigen/Dense>

class Earth
{

public:
    static constexpr double EARTH_ROTATION_RATE = 7.2921150e-5; // rad/s

    static constexpr double MU = 3.986004418e14; // m^3 / s^2

    static constexpr double J2 = 1.75553e19; // m^5 / s^2

    static constexpr double EARTH_EQUATOR_R = 6378137; // m

    static constexpr double EARTH_POLAR_R = 6356752.314; // m

    struct Geodetic 
    {
        double latitude;
        double longitude;
        double altitude;
    };

    static Eigen::Vector3d get_fictional_forces(const Eigen::Vector3d& ecef_position, 
        const Eigen::Vector3d& ecef_velocity);

    static Eigen::Vector3d get_gravity(const Eigen::Vector3d& ecef_position);

    static Eigen::Vector3d get_J2_gravity(const Eigen::Vector3d& ecef_position);

    static double get_altitude(const Eigen::Vector3d& ecef_position);

    static Geodetic ECEF2LLA(const Eigen::Vector3d& ecef_position);

    static Eigen::Vector3d LLA2ECEF(const Geodetic& LLA);

    static Eigen::Matrix3d getENU(double longitude, double latitude);

    static Eigen::Matrix3d getNED(double longitude, double latitude);

    static Eigen::Vector3d getUp(const Eigen::Vector3d& ecef);

    static Eigen::Matrix3d getENU(const Eigen::Vector3d& ecef);

    static void getGeodeticFrame(const Eigen::Vector3d& ecef_position, Geodetic& lla, double* ENU);

    static std::array<double,6> ECI2Keplerian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity);
};

#endif