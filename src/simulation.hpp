#ifndef SIMULATION_H
#define SIMULATION_H

#include "dynamic_array.hpp"
#include "ode.hpp"
#include "vehicle.hpp"
#include "track_dynamics.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <array>
#include <string>
#include <vector>


namespace Simulation
{
    struct DerivedQuantities
    {
        double g_force;
        double flight_angle;
        double mach;
        double altitude;
    };

    struct SimulationResult
    {
        std::vector<double> times;
        std::vector<State> ecef_states;
        std::vector<DerivedQuantities> derived;
        std::vector<dynamic_array<double>> other_data;
        double start_ramjet;
        double start_rocket;
        std::array<double, 6> final_orbital_elements;
    };

    std::vector<DerivedQuantities> get_derived(const std::vector<double>& times, 
        const std::vector<State>& ecef, const Atmosphere& atm);

    json load_config(const std::string& fileName);

    SimulationResult run(json& config);

}


#endif