#include "earth.hpp"

#include <array>
#include <cmath>

Eigen::Vector3d Earth::get_fictional_forces(const Eigen::Vector3d& ecef_position, 
        const Eigen::Vector3d& ecef_velocity)
{
    double x = -Earth::EARTH_ROTATION_RATE*(Earth::EARTH_ROTATION_RATE*ecef_position[0] - 2.0*ecef_velocity[1]);
    double y = -Earth::EARTH_ROTATION_RATE*(Earth::EARTH_ROTATION_RATE*ecef_position[1] + 2.0*ecef_velocity[0]);

    return Eigen::Vector3d(x,y,0);
}

Eigen::Vector3d Earth::get_gravity(const Eigen::Vector3d& ecef_position)
{
    double r2 = ecef_position.dot(ecef_position);
    r2 = -Earth::MU/(r2*sqrt(r2));
    return ecef_position*r2;
}

Eigen::Vector3d Earth::get_J2_gravity(const Eigen::Vector3d& ecef_position)
{
    double x2 = ecef_position.x()*ecef_position.x();
    double y2 = ecef_position.y()*ecef_position.y();
    double z2 = ecef_position.z()*ecef_position.z();
    double w2 = x2 + y2;
    double inv_r2 = 1.0/(w2 + z2);

    double numxy = (6.0*z2 + 1.5*w2)*inv_r2;
    double numz = numxy - 3.0;

    double inv_r3 = sqrt(inv_r2)*inv_r2;
    double J2_scaled = Earth::J2*inv_r3*inv_r2;
    double MU_scaled = Earth::MU*inv_r3;

    Eigen::Vector3d g;

    g(0) = -ecef_position.x()*(MU_scaled + J2_scaled*numxy);
    g(1) = -ecef_position.y()*(MU_scaled + J2_scaled*numxy);
    g(2) = -ecef_position.z()*(MU_scaled + J2_scaled*numz);

    return g;
}

constexpr double e2 = 6.6943799901377997e-3;  //WGS-84 first eccentricity squared
constexpr double a1 = 4.2697672707157535e+4;  //a1 = a*e2
constexpr double a2 = 1.8230912546075455e+9;  //a2 = a1*a1
constexpr double a3 = 1.4291722289812413e+2;  //a3 = a1*e2/2
constexpr double a4 = 4.5577281365188637e+9;  //a4 = 2.5*a2
constexpr double a5 = 4.2840589930055659e+4;  //a5 = a1+a3
constexpr double a6 = 9.9330562000986220e-1;  //a6 = 1-e2

double Earth::get_altitude(const Eigen::Vector3d& ecef_position)
{
    double zp = fabs(ecef_position.z());
    double z2 = zp * zp;
    double w2 = ecef_position.x()*ecef_position.x() + ecef_position.y()*ecef_position.y();
    double w = sqrt(w2);
    double r2 = 1.0 / (w2 + z2);
    double r_inv = sqrt(r2);
    double s2 = z2 * r2;
    double c2 = w2 * r2;
    double u = a2 * r_inv;
    double v = a3 - a4 * r_inv;
    double s, c, ss;
    if (c2 > 0.3) {
        s = (zp * r_inv) * (1.0 + c2 * (a1 + u + s2 * v) * r_inv);
        ss = s * s;
        c = sqrt(1.0 - ss);
    }
    else {
        c = (w * r_inv) * (1.0 - s2 * (a5 - u - c2 * v) * r_inv);
        ss = 1.0 - c * c;
        s = sqrt(ss);
    }
    double g = 1.0 - e2 * ss;
    double rg = Earth::EARTH_EQUATOR_R / sqrt(g);
    double rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    double f = c * u + s * v;
    double m = c * v - s * u;
    double p = m / (rf / g + f);
    return f + m * p * 0.5;
}

Earth::Geodetic Earth::ECEF2LLA(const Eigen::Vector3d& ecef_position)
{
    Earth::Geodetic lla;
    lla.longitude = atan2(ecef_position.y(), ecef_position.x());
    double zp = fabs(ecef_position.z());
    double z2 = zp * zp;
    double w2 = ecef_position.x()*ecef_position.x() + ecef_position.y()*ecef_position.y();
    double w = sqrt(w2);
    double r2 = 1.0 / (w2 + z2);
    double r_inv = sqrt(r2);
    double s2 = z2 * r2;
    double c2 = w2 * r2;
    double u = a2 * r_inv;
    double v = a3 - a4 * r_inv;
    double s, c, ss;
    if (c2 > 0.3) {
        s = (zp * r_inv) * (1.0 + c2 * (a1 + u + s2 * v) * r_inv);
        lla.latitude = asin(s);
        ss = s * s;
        c = sqrt(1.0 - ss);
    }
    else {
        c = (w * r_inv) * (1.0 - s2 * (a5 - u - c2 * v) * r_inv);
        lla.latitude = acos(c);
        ss = 1.0 - c * c;
        s = sqrt(ss);
    }
    double g = 1.0 - e2 * ss;
    double rg = Earth::EARTH_EQUATOR_R / sqrt(g);
    double rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    double f = c * u + s * v;
    double m = c * v - s * u;
    double p = m / (rf / g + f);
    lla.latitude += p;
    lla.altitude = f + m * p * 0.5;
    if( ecef_position.z() < 0.0 ) {
        lla.latitude *= -1.0;     //Lat
    }
    return lla;
}

void Earth::getGeodeticFrame(const Eigen::Vector3d& ecef_position, Earth::Geodetic& lla, double* ENU)
{
    lla.longitude = atan2(ecef_position.y(), ecef_position.x());
    const auto zp = fabs(ecef_position.z());
    const auto z2 = zp * zp;
    const auto w2 = ecef_position.x()*ecef_position.x() + ecef_position.y()*ecef_position.y();
    const auto w = sqrt(w2);
    const auto r2 = 1.0 / (w2 + z2);
    const auto r_inv = sqrt(r2);
    const auto s2 = z2 * r2;
    const auto c2 = w2 * r2;
    auto u = a2 * r_inv;
    auto v = a3 - a4 * r_inv;
    double s, c, ss;
    if (c2 > 0.3) {
        s = (zp * r_inv) * (1.0 + c2 * (a1 + u + s2 * v) * r_inv);
        lla.latitude = asin(s);
        ss = s * s;
        c = sqrt(1.0 - ss);
    }
    else {
        c = (w * r_inv) * (1.0 - s2 * (a5 - u - c2 * v) * r_inv);
        lla.latitude = acos(c);
        ss = 1.0 - c * c;
        s = sqrt(ss);
    }
    double g = 1.0 - e2*ss;
    double rg = Earth::EARTH_EQUATOR_R / sqrt(g);
    double rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    double f = c * u + s * v;
    double m = c * v - s * u;
    double p = m / (rf / g + f);
    lla.latitude += p;
    lla.altitude = f + m * p * 0.5;
    if( ecef_position.z() < 0.0 ) {
        lla.latitude *= -1.0;     //Lat
    }
    
    // 
    ENU[0] = -ecef_position.y()/w;
    ENU[1] = ecef_position.x()/w;
    // 
    ENU[3] = -ENU[0]*s;
    ENU[4] = ENU[1]*s;
    ENU[5] = c;
    ENU[6] = ENU[1]*c;
    ENU[7] = -ENU[0]*c;
    ENU[8] = s;
}

Eigen::Vector3d Earth::LLA2ECEF(const Geodetic& LLA)
{
    double s = sin(LLA.latitude);
    auto n = EARTH_EQUATOR_R/sqrt(1.0 - e2*s*s);
    auto tmp = (n + LLA.altitude)*sqrt(1 - s*s);
    Eigen::Vector3d ecef;
    ecef(0) = tmp*cos(LLA.longitude);
    ecef(1) = tmp*sin(LLA.longitude);
    ecef(2) = s*(n*(1.0 - e2) + LLA.altitude);
    return ecef;
}

Eigen::Matrix3d Earth::getENU(double longitude, double latitude)
{
    double east_x = -sin(longitude);
    double east_y = cos(longitude);

    double up_z = sin(latitude);
    double north_z = sqrt(1.0 - up_z*up_z);
    double up_x = east_y*north_z;
    double up_y = -east_x*north_z;
    double north_x = east_x*up_z;
    double north_y = -east_y*up_z;
    
    Eigen::Matrix3d ENU;
    ENU << east_x, north_x, up_x, east_y, north_y, up_y, 0.0, north_z, up_z;
    return ENU;
}

Eigen::Matrix3d Earth::getNED(double longitude, double latitude)
{
    double east_x = -sin(longitude);
    double east_y = cos(longitude);

    double up_z = sin(latitude);
    double north_z = sqrt(1.0 - up_z*up_z);
    double up_x = east_y*north_z;
    double up_y = -east_x*north_z;
    double north_x = east_x*up_z;
    double north_y = -east_y*up_z;
    
    Eigen::Matrix3d NED;
    NED << north_x, east_x, -up_x, north_y, east_y, -up_y, north_z, 0.0, -up_z;
    return NED;
}

#include <iostream>
/*
Eigen::Vector3d Earth::getUp(const Eigen::Vector3d& ecef)
{
    const double z2 = ecef.z()*ecef.z();
    const double p2 = ecef.x()*ecef.x() + ecef.y()*ecef.y();

    const double inv_p2 = 1.0/p2;

    constexpr double KAPPA = EARTH_EQUATOR_R*EARTH_EQUATOR_R - EARTH_POLAR_R*EARTH_POLAR_R;;
    constexpr double KAPPA_SQ = KAPPA*KAPPA;
    double k = KAPPA_SQ*inv_p2;
    double beta_sq = z2*inv_p2;

    double z_prime_sq = z2/(z2 + p2);
    double p_prime_sq = 1.0 - z_prime_sq;
    std::cout << z_prime_sq << std::endl;
    for(int i = 0; i < 1; i++)
    {
        double alpha = sqrt(k*p_prime_sq/(1.0 - e2*z_prime_sq));
        double nu = 1.0 - alpha;
        double D = nu*(beta_sq + nu*nu);
        double E = nu*beta_sq;
        double num = D*z_prime_sq - E;
        double den = D + E*a6*alpha;
        z_prime_sq -= num/den;
        p_prime_sq = 1.0 - z_prime_sq;
        std::cout << z_prime_sq << std::endl;
    } 

    double p_scale = sqrt(p_prime_sq*inv_p2);
    double z_prime = sqrt(z_prime_sq);
    
    return Eigen::Vector3d(ecef.x()*p_scale, ecef.y()*p_scale, z_prime);
}
*/

Eigen::Vector3d Earth::getUp(const Eigen::Vector3d& ecef)
{
    const double z2 = ecef.z()*ecef.z();
    const double p2 = ecef.x()*ecef.x() + ecef.y()*ecef.y();

    const double inv_p2 = 1.0/p2;

    constexpr double KAPPA_SQ = a1*a1;
    const double k = KAPPA_SQ*inv_p2;
    const double beta_sq = z2*inv_p2;

    double z_prime_sq = z2/(z2 + p2);
    double p_prime_sq = 1.0 - z_prime_sq;
    for(int i = 0; i < 4; i++)
    {
        double nu = 1.0 - sqrt(k*p_prime_sq/(1.0 - e2*z_prime_sq));
        z_prime_sq = beta_sq*p_prime_sq/(nu*nu);
        p_prime_sq = 1.0 - z_prime_sq;
    } 

    double p_scale = sqrt(p_prime_sq*inv_p2);
    double z_prime = sqrt(z_prime_sq);
    
    return Eigen::Vector3d(ecef.x()*p_scale, ecef.y()*p_scale, z_prime);
}
// Look up table

Eigen::Matrix3d Earth::getENU(const Eigen::Vector3d& ecef)
{
    Eigen::Matrix3d ENU;

    const double z2 = ecef.z()*ecef.z();
    const double p2 = ecef.x()*ecef.x() + ecef.y()*ecef.y();

    const double inv_p2 = 1.0/p2;
    const double inv_r2 = 1.0/(z2 + p2);
    const double inv_p = sqrt(inv_p2);

    constexpr double KAPPA_SQ = a1*a1;
    const double k = KAPPA_SQ*inv_p2;
    const double beta_sq = z2*inv_p2;

    double z_prime_sq = z2*inv_r2;
    double p_prime_sq = 1.0 - z_prime_sq;
    for(int i = 0; i < 4; i++)
    {
        double nu = 1.0 - sqrt(k*p_prime_sq/(1.0 - e2*z_prime_sq));
        z_prime_sq = beta_sq*p_prime_sq/(nu*nu);
        p_prime_sq = 1.0 - z_prime_sq;
    } 

    double p_scale = sqrt(p_prime_sq*inv_p2);
    double z_prime = sqrt(z_prime_sq);

    ENU(0,0) = -ecef.y()*inv_p;
    ENU(1,0) = ecef.x()*inv_p;
    ENU(2,0) = 0.0;

    ENU(0,2) = ecef.x()*p_scale;
    ENU(1,2) = ecef.y()*p_scale;
    ENU(2,2) = z_prime;

    ENU(0,1) = ENU(1,2)*ENU(2,0) - ENU(2,2)*ENU(1,0);
    ENU(1,1) = ENU(2,2)*ENU(0,0) - ENU(0,2)*ENU(2,0);
    ENU(2,1) = ENU(0,2)*ENU(1,0) - ENU(1,2)*ENU(0,0);

    return ENU;
}

std::array<double,6> Earth::ECI2Keplerian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity)
{
    constexpr double TWOPI = 6.2831853071795864769;
    std::array<double, 6> oe{};
    Eigen::Vector3d h = position.cross(velocity);

    // std::array<double,2> n = {h[1],-h[0]}; // z is implicit 0
    double v2 = velocity.dot(velocity);
    double r_inv = 1.0 / position.norm();
    double rv = position.dot(velocity);
    const double invMU = 1.0 / MU;
    double tmp1 = v2*invMU - r_inv;
    double tmp2 = rv*invMU;
    Eigen::Vector3d e  = position*tmp1 + velocity*tmp2;
    double egy = v2*0.5 - MU*r_inv;

    oe[0] = -MU / (2*egy);
    oe[1] = e.norm();
    double inv_e = 1 / oe[1];
    double nmag = h[0]*h[0] + h[1]*h[1];
    oe[2] = acos(h[2] / sqrt(nmag + h[2]*h[2]));

    if (fabs(oe[2]) > 1e-9) {
        nmag = 1.0 / sqrt(nmag);
        oe[3] = acos(h[1]*nmag);
        if (h[0] > 0) 
        {
            oe[3] = TWOPI - oe[3];
        }

        oe[4] = acos((h[1]*e[0] - h[0]*e[1])*nmag*inv_e);
        if (e[2] < 0) 
        {
            oe[4] = TWOPI - oe[4];
        }
    }

    oe[5] = acos(e.dot(position)*r_inv*inv_e);
    if (rv < 0) 
    {
        oe[5] = TWOPI - oe[5];
    }
    return oe;
}