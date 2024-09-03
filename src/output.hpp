#ifndef OUTPUT_H
#define OUTPUT_H

#include "ode.hpp"
#include "simulation.hpp"
#include "track_dynamics.hpp"

#include <fstream>
#include <string>
#include <mutex>

class Logging 
{
    static std::ofstream _file;

    static const std::array<std::string,3> TYPE_LABELS;

    static std::mutex _lock;

public:

    enum LOG_TYPE
    {
        INFO = 0,
        WARN,
        ERROR
    };

    static void init(const std::string& filename = "main.log")
    {
        if(!_file.is_open())
        {
            _file.open(filename, std::ofstream::trunc);
        }
    }

    static void log(LOG_TYPE type, int level, const std::string& message);

    static void debug_log(int level, const std::string& message);
    
};

namespace output
{

    void write_track_csv(const std::string& filename, const Track& track);

    template<unsigned NSTATES>
    void write_recording_csv(const std::string& filename, const Recording<NSTATES>& recording);

    void write_simulation_result(const std::string& filename, const Simulation::SimulationResult& result);
}

#endif