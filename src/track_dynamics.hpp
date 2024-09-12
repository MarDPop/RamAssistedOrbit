#ifndef TRACK_DYNAMICS_H
#define TRACK_DYNAMICS_H

#include "constants.hpp"
#include "state_6dof.hpp"

#include <vector>

struct Track
{
    std::vector<std::array<double,4>> x;

    std::vector<double> t;
};

namespace track_dynamics
{
    Track generate_track(double max_g, double friction_coef, double linear_friction_coef,
        double launch_exit_speed, double launch_exit_angle, double dt, double dt_record);

    State_6DOF generate_exit_ecef_state(double launch_longitude, 
        double launch_latitude, double launch_altitude, double launch_heading, double launch_speed, double launch_angle);

    std::vector<State_6DOF> convert_track_states(const Track& track, double mass, double launch_longitude, 
        double launch_latitude, double launch_altitude, double launch_heading, double launch_speed, double launch_angle);

}

#endif