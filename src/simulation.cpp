#include "simulation.hpp"

#include "atmosphere.hpp"
#include "constants.hpp"
#include "earth.hpp"
#include "ode.hpp"
#include "output.hpp"
#include "track_dynamics.hpp"
#include "vehicle.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#define DEG2RAD 0.01745329251994

namespace Simulation
{

json load_config(const std::string& fileName)
{
    json config;
    std::ifstream file(fileName);

    if (!file.is_open()) {
        throw std::runtime_error("cannot find configuration file");
    }

    file >> config;
    file.close();
    return config;
}

std::vector<DerivedQuantities> get_derived(const std::vector<double>& times, 
    const std::vector<State_6DOF>& states, const Atmosphere& atm)
{
    std::vector<DerivedQuantities> derived(times.size());
    for(auto i = 0u; i < times.size(); i++)
    {
        const auto& state = states[i].body;
        double speed = state.velocity.norm();
        const auto lla = Earth::ECEF2LLA(state.position);
        derived[i].altitude = lla.altitude;
        auto ENU = Earth::getENU(lla.longitude, lla.latitude);
        auto vertical_speed = ENU.col(2).dot(state.velocity);
        derived[i].flight_angle = asin(vertical_speed / speed);
        Air air;
        atm.set_air(std::max(0.0,lla.altitude), air);
        derived[i].mach = speed * air.inv_sound_speed;
        if(i > 0 && i < (times.size() - 1))
        {
            Eigen::Vector3d acceleration = (states[i+1].body.velocity - states[i-1].body.velocity)/(GForce::G*(times[i+1] - times[i-1]));
            acceleration += (state.position / state.position.norm()); 
            derived[i].g_force = acceleration.norm();
        }
        else
        {
            derived[i].g_force = 1.0;
        }
    }
    
    return derived;
}

std::array<double,6> get_final_orbital_elements(const State_6DOF& state)
{
    Eigen::Vector3d velocity_eci = state.body.velocity;
    velocity_eci(0) -= Earth::EARTH_ROTATION_RATE*state.body.position(1);
    velocity_eci(1) += Earth::EARTH_ROTATION_RATE*state.body.position(0);
    return Earth::ECI2Keplerian(state.body.position, velocity_eci);
}

SimulationResult run(json& config)
{
    std::cout << "**** Running Track States **** \n";
    SimulationResult result;

    constexpr double dt = 0.001;
    constexpr double dt_record = 0.1;

    auto track = track_dynamics::generate_track(config["TRACK"]["MAX_G"].template get<double>(), 
        config["TRACK"]["FRICTION"].template get<double>(), 
        config["TRACK"]["LINEAR_FRICTION"].template get<double>(), 
        config["TRACK"]["EXIT_SPEED"].template get<double>(), 
        config["TRACK"]["EXIT_ANGLE"].template get<double>()*DEG2RAD, 
        dt, dt_record);

    output::write_track_csv(config["OUTPUT_DIRECTORY"].template get<std::string>() + "track.dat", track);

    auto track_end_state = track_dynamics::generate_exit_ecef_state(config["TRACK"]["LAUNCH_LONGITUDE"].template get<double>()*DEG2RAD, 
        config["TRACK"]["LAUNCH_LATITUDE"].template get<double>()*DEG2RAD, 
        config["TRACK"]["LAUNCH_ALTITUDE"].template get<double>(), 
        config["TRACK"]["LAUNCH_HEADING"].template get<double>()*DEG2RAD, 
        config["TRACK"]["EXIT_SPEED"], config["TRACK"]["EXIT_ANGLE"].template get<double>()*DEG2RAD);

    std::cout << "**** Track End State **** \n";
    std::cout << track_end_state.body.position << "\n";
    std::cout << track_end_state.body.velocity << "\n";
    std::cout << "**** Setup Vehicle Simulation ****" << std::endl;

    const double startMass = config["INERTIAL_PROPERTIES"]["MASS_0"].template get<double>();

    result.ecef_states = track_dynamics::convert_track_states(track, startMass, 
        config["TRACK"]["LAUNCH_LONGITUDE"].template get<double>()*DEG2RAD, 
        config["TRACK"]["LAUNCH_LATITUDE"].template get<double>()*DEG2RAD, 
        config["TRACK"]["LAUNCH_ALTITUDE"].template get<double>(), 
        config["TRACK"]["LAUNCH_HEADING"].template get<double>()*DEG2RAD, 
        config["TRACK"]["EXIT_SPEED"].template get<double>(), 
        config["TRACK"]["EXIT_ANGLE"].template get<double>()*DEG2RAD);

    result.times = track.t;
    result.start_ramjet = track.t.back();

    for(unsigned i = 0; i < track.t.size(); i++) {
        result.other_data.emplace_back(0);
    }

    // Ramjet Stage
    AtmosphereLinearTable atm = AtmosphereLinearTable::create(AtmosphereLinearTable::STD_ATMOSPHERES::US_1976, 100);

    constexpr double min_empty = 0.01;
    InertialProperties I(startMass*min_empty, 
        {config["INERTIAL_PROPERTIES"]["IXX_0"].template get<double>()*min_empty,
        config["INERTIAL_PROPERTIES"]["IYY_0"].template get<double>()*min_empty,
        config["INERTIAL_PROPERTIES"]["IZZ_0"].template get<double>()*min_empty}, Eigen::Vector3d(0,0,0),
        startMass, 
        {config["INERTIAL_PROPERTIES"]["IXX_0"].template get<double>(),
        config["INERTIAL_PROPERTIES"]["IYY_0"].template get<double>(),
        config["INERTIAL_PROPERTIES"]["IZZ_0"].template get<double>()}, Eigen::Vector3d(0,0,0));

    AerodynamicBasicCoefficients::Coef coef;
    coef.CD0 = config["RAMJET"]["AERODYNAMICS"]["CD0"].template get<double>();
    coef.K = config["RAMJET"]["AERODYNAMICS"]["INDUCED_DRAG_K"].template get<double>();
    coef.alpha0 = config["RAMJET"]["AERODYNAMICS"]["ALPHA0"].template get<double>();
    coef.CL_alpha = config["RAMJET"]["AERODYNAMICS"]["CL_ALPHA"].template get<double>();
    coef.CM_alpha = config["RAMJET"]["AERODYNAMICS"]["CM_ALPHA"].template get<double>();
    coef.CN_beta = config["RAMJET"]["AERODYNAMICS"]["CN_BETA"].template get<double>();
    coef.stall_angle = config["RAMJET"]["AERODYNAMICS"]["STALL_ANGLE"].template get<double>();
    coef.dCD_dEl = config["RAMJET"]["AERODYNAMICS"]["DCD_DEL"].template get<double>();
    coef.dCM_dEl = config["RAMJET"]["AERODYNAMICS"]["DCM_DEL"].template get<double>();
    coef.max_elevator_deflection = config["RAMJET"]["AERODYNAMICS"]["MAX_EL_DEFLECTION"].template get<double>();
    coef.A_ref = config["RAMJET"]["AERODYNAMICS"]["AREF"].template get<double>();
    coef.L_ref = config["RAMJET"]["AERODYNAMICS"]["LREF"].template get<double>();

    const double cruise_lift2drag = 4*(config["RAMJET"]["CRUISE"]["MACH"].template get<double>() + 3)/
         config["RAMJET"]["CRUISE"]["MACH"].template get<double>();

    auto ramjet = RamjetVariableInlet::create(config["RAMJET"]["THRUST2WEIGHT"].template get<double>(), 
        config["RAMJET"]["CRUISE"]["ALTITUDE"].template get<double>(), 
        config["RAMJET"]["CRUISE"]["MACH"].template get<double>(), 
        config["RAMJET"]["CRUISE"]["MACH"].template get<double>(), cruise_lift2drag, 
        startMass, 
        config["RAMJET"]["THRUST_MARGIN"].template get<double>());

    ramjet.stop();

    double thrust2weight = config["RAMJET"]["THRUST2WEIGHT"].template get<double>();
    double thrust = startMass*thrust2weight*GForce::G;
    auto ramjet2 = RamjetFixedInlet::create(config["RAMJET"]["CRUISE"]["MACH"].template get<double>(),
        config["RAMJET"]["CRUISE"]["ALTITUDE"].template get<double>(), 
        thrust,
        8.0);
    
    RamjetVehicle rVehicle(I, atm, ramjet2, coef);

    rVehicle.set_control_values(config["RAMJET"]["CONTROL"]["K1"].template get<double>(), 
        config["RAMJET"]["CONTROL"]["C1"].template get<double>(), 
        config["RAMJET"]["CONTROL"]["MAX_ALPHA"].template get<double>(), 
        config["RAMJET"]["CONTROL"]["MIN_ALPHA"].template get<double>(), 
        config["RAMJET"]["CONTROL"]["ALPHA_K"].template get<double>(), 
        config["RAMJET"]["CRUISE"]["ALTITUDE"].template get<double>(),
        config["RAMJET"]["CRUISE"]["MACH"].template get<double>());

    //ODE_HUEN_EULER<RamjetVehicle> ode(rVehicle);
    int odeType = config["ODE"]["TYPE"].template get<int>();
    ODE<RamjetVehicle>* ode;
    switch(odeType)
    {
        case 1:
            ode = new ODE_HUEN_EULER<RamjetVehicle>(rVehicle);
            break;
        case 2:
            ode = new ODE_RK4<RamjetVehicle>(rVehicle);
            break;
        case 3: 
            ode = new ODE_RK45<RamjetVehicle>(rVehicle);
            break;
        default:
            ode = new ODE<RamjetVehicle>(rVehicle);
    }

    RunOptions<14> options;

    options.initial_state = track_end_state.x;
    options.initial_state[13] = startMass;
    options.initial_time = 0;
    options.final_time = config["ODE"]["MAX_RUNTIME"].template get<double>();
    options.recording_interval = config["ODE"]["RECORDING_INTERVAL"].template get<double>();
    options.max_steps = 10000000;
    options.timestep.max_stepsize = config["ODE"]["MAX_STEP_SIZE"].template get<double>();
    options.timestep.min_stepsize = config["ODE"]["MIN_STEP_SIZE"].template get<double>();
    options.timestep.inv_absolute_error = std::array<double, 14> { 
        1e3, 1e3, 1e3, 1e2, 1e2, 1e2, 1e5, 1e5, 1e5, 1e5, 1e3, 1e3, 1e3, 1e5
    };

    std::cout << "**** Running Vehicle Simulation ****" << std::endl;
    try 
    {
        ode->run(options);
    } 
    catch (std::exception& e) 
    {
        std::cerr << e.what() << std::endl;
    } 
    catch(...) 
    {
        std::cout << "unknown error occured" << std::endl;
        std::cout << "num entries: " << ode->recording().times.size() << std::endl;
    }
    std::cout << "**** Simulation Finished ****" << std::endl;

    if(ode->dynamics_stopped())
    {
        std::cout << "ODE simulation stopped prematurely" << std::endl;
        std::cout << "Dynamics invalid" << std::endl;
    } 
    else if(ode->options_stopped())
    {
        std::cout << "ODE simulation stopped correctly" << std::endl;
    } 
    else 
    {
        std::cout << "ODE simulation ran out of time" << std::endl;
    }

    const unsigned nData = rVehicle.num_data();
    for(auto i = 0u; i < ode->recording().times.size(); i++)
    {
        result.times.push_back(ode->recording().times[i] + result.start_ramjet);
        State_6DOF state;
        state.x = ode->recording().states[i];
        result.ecef_states.push_back(state);
        dynamic_array<double> other(nData);
        for(auto j = 0u; j < nData; j++) {
            other[j] = ode->recording().output[i][j];
        }
        result.other_data.push_back(other);
    }

    std::cout << "entries: " << ode->recording().times.size() << std::endl;
    
    // Finish
    result.final_orbital_elements = get_final_orbital_elements(result.ecef_states.back());

    std::cout << "Orbital El: ";
    for(auto e : result.final_orbital_elements)
    {
        std::cout << e << " ";
    }
    std::cout << std::endl;

    result.derived = get_derived(result.times, result.ecef_states, atm);

    return result;
}

}// end namespace