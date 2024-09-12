#ifndef VEHICLE_ODE_HPP
#define VEHICLE_ODE_HPP

#include "vehicle.hpp"

struct Vehicle_Recording 
{
    std::vector<State_6DOF> states;
    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> moments;
    std::vector<std::vector<double>> other_values;
    std::vector<double> times;
};

class Vehicle_ODE 
{
public:
    
    struct Vehicle_ODE_Options
    {
        double max_position_error;
        double max_velocity_error;
        double max_angle_error;
        double max_angular_velocity_error;
        double min_time_step;
        double max_time_step;
        double recording_interval;
        double fake_rotation_damping;
    };

    Vehicle_ODE_Options options;

    bool (*stop)(const Vehicle&) =  [](const Vehicle&){ return false; };

    Vehicle_ODE(Vehicle& v) : _vehicle(v) {}

    const Vehicle_Recording& get_recording() const { return _recording; }

    void reset_recording() { _recording = Vehicle_Recording(); }

    void reset_time(double time = 0.0) { _t = time; }

    void run(double run_time)
    {
        double t_record = _t;
        _dt = options.min_time_step;
        const double t_end = _t + run_time;
        while (_t < t_end && !stop(_vehicle))
        {
            if(_t >= t_record)
            {
                _recording.times.push_back(_t);
                _recording.states.push_back(_vehicle.get_state());
                _recording.forces.push_back(_vehicle.get_body_force());
                _recording.moments.push_back(_vehicle.get_body_moment());
                _recording.other_values.push_back(_vehicle.get_other_values());
                t_record = _recording.times.size()*options.recording_interval;
            }
            update();
            _t += _dt;
        }
    }



    virtual void update() 
    {
        _vehicle.update(_t);
        // add fake damping for better stability
        _vehicle._body_moment -= (_vehicle._state.body.angular_velocity*options.fake_rotation_damping); 
        // remember angular_velocity and moment is in body frame

        for(int iter = 0; iter < 4; ++iter)
        {

            _vehicle.update(_t);
        }
    }

protected:

    Vehicle& _vehicle; 

    Vehicle_Recording _recording;

    double _t;

    double _dt;

};



#endif