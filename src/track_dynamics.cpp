#include "track_dynamics.hpp"

#include "constants.hpp"
#include "earth.hpp"
#include "functions.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

Track track_dynamics::generate_track(double max_g, double friction_coef, double linear_friction_coef,
        double launch_exit_speed, double launch_exit_angle, double dt, double dt_record)
{
    double position[2] = {0.0,0.0};
    double velocity[2] = {cos(launch_exit_angle)*launch_exit_speed, 
        sin(launch_exit_angle)*launch_exit_speed};

    double acceleration[2];

    double direction[2];

    const double maxA = max_g*GForce::G;
    const double time_to_coast = 2.0;
    const double dt2 = dt*0.5;

    double time = 0;
    double time_record = -dt_record;
    bool pullup = true;
    bool coast = false;
    double t_coast = 0;

    Track track;
    track.t.push_back(time);
    track.x.push_back({position[0], position[1], velocity[0], velocity[1]});

    for(int iter = 0; iter < 1000000; iter++)
    {
        if(time < time_record)
        {
            track.t.push_back(time);
            track.x.push_back({position[0], position[1], velocity[0], velocity[1]});
            time_record -= dt_record;
        }

        if(velocity[0] < 0)
        {
            track.t.push_back(time);
            track.x.push_back({position[0], position[1], velocity[0], velocity[1]});
            break;
        }

        double speed = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1]);

        direction[0] = velocity[0] / speed;
        direction[1] = velocity[1] / speed;

        if(pullup)
        {
            double a_up = (max_g - direction[0])*GForce::G;
            double drag = linear_friction_coef*speed + friction_coef*max_g;
            double a_forward = -drag - direction[1]*GForce::G;
            
            acceleration[0] = a_forward*direction[0] - a_up*direction[1];
            acceleration[1] = a_forward*direction[1] + a_up*direction[0];

            position[0] -= dt*(velocity[0] - acceleration[0]*dt2);
            position[1] -= dt*(velocity[1] - acceleration[1]*dt2);

            velocity[0] -= dt*acceleration[0];
            velocity[1] -= dt*acceleration[1];

            pullup = direction[1] > 0;
            if(!pullup)
            {
                t_coast = time - time_to_coast;
                coast = true;
            }
        }
        else
        {
            double a = maxA;
            if(coast)
            {
                a = -(linear_friction_coef*speed + friction_coef);
                if(time < t_coast)
                {
                    coast = false;
                }
            }
            position[0] -= dt*(velocity[0] - a*dt2);
            velocity[0] -= dt*a;
        }

        time -= dt;       
    }

    for(auto i = 0u; i < track.t.size(); i++)
    {
        track.t[i] -= time;
        track.x[i][0] -= position[0];
        track.x[i][1] -= position[1];
    }
    std::reverse(track.t.begin(), track.t.end());
    std::reverse(track.x.begin(), track.x.end());

    return track;
}

State_6DOF track_dynamics::generate_exit_ecef_state(double launch_longitude, 
        double launch_latitude, double launch_altitude, double launch_heading, double launch_speed, double launch_angle)
{
    State_6DOF ecef_state;
    ecef_state.x.fill(0.0);
    // get ECEF position
    Earth::Geodetic lla = {launch_longitude, launch_latitude, launch_altitude};
    ecef_state.body.position = Earth::LLA2ECEF(lla);

    // Get rotation of the heading and track exit angle
    Eigen::Matrix3d ENU_rotation; // rotation of the ENU frame
    functions::ZY_rotation(launch_heading, -launch_angle, ENU_rotation.data()); // col major

    // convert to vehicle principal axis (z down)
    Eigen::Matrix3d vehicle_ENU_orientation;
    vehicle_ENU_orientation << ENU_rotation(0,0), -ENU_rotation(0,1), -ENU_rotation(0,2),
        ENU_rotation(1,0), -ENU_rotation(1,1), -ENU_rotation(1,2),
        ENU_rotation(2,0), -ENU_rotation(2,1), -ENU_rotation(2,2);

    // Get local tangent plane orientation
    Eigen::Matrix3d ENU2ECEF = Earth::getENU(launch_longitude, launch_latitude);

    // multiply column major orientation to convert to ECEF frame then transpose to be row major
    Eigen::Matrix3d vehicle_ecef_orientation = (ENU2ECEF*vehicle_ENU_orientation).transpose(); // row major

    // convert to quaternion
    ecef_state.body.orientation = Eigen::Quaterniond(vehicle_ecef_orientation);
    
    ecef_state.body.velocity = (vehicle_ecef_orientation.row(0))*launch_speed;

    return ecef_state;
}

std::vector<State_6DOF> track_dynamics::convert_track_states(const Track& track, double mass, double launch_longitude, 
        double launch_latitude, double launch_altitude, double launch_heading, double launch_speed, double launch_angle)
{
    std::vector<State_6DOF> states;

    State_6DOF ecef_end_state = track_dynamics::generate_exit_ecef_state(launch_longitude, launch_latitude, launch_altitude, 
        launch_heading, launch_speed, launch_angle);

    Eigen::Matrix3d NED = Earth::getNED(launch_longitude, launch_latitude);

    Eigen::Vector3d east = NED.col(1);
    Eigen::Vector3d north = NED.col(0);
    Eigen::Vector3d track_x_axis = east*cos(launch_heading) + north*sin(launch_heading);
    Eigen::Vector3d track_z_axis = NED.col(2);

    double x_final = track.x.back()[0];
    double z_final = track.x.back()[1];

    // assume flat LTP ( ~ 10 meter drop in 8km track )
    double z_old = 0;
    for(const auto& x : track.x)
    {
        State_6DOF state;
        state.body.position = ecef_end_state.body.position + (x[0] - x_final)*track_x_axis - (x[1] - z_final)*track_z_axis;
        state.body.velocity = track_x_axis*x[2] - track_z_axis*x[3];
        double speed = sqrt(x[2]*x[2] + x[3]*x[3]);
        Eigen::Vector3d forward = state.body.velocity / speed;
        Eigen::Vector3d down = (track_z_axis*x[2] + track_x_axis*x[3]) / speed;
        Eigen::Vector3d side = down.cross(forward);
        Eigen::Matrix3d rotm;
        // orientation is x forward, y to right, z down
        rotm << forward(0), forward(1), forward(2), side(0), side(1), side(2), down(0), down(1), down(2); 
        state.body.orientation = Eigen::Quaterniond(rotm);
        if(states.size() > 0)
        {
            double omega = x[3] - z_old;
            state.body.angular_velocity = side*omega;
            z_old = x[3];
        }
        else
        {
            state.body.angular_velocity.setZero();
        }
        
        state.body.mass = mass;

        states.push_back(state);
    }

    return states;
}