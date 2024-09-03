#include "output.hpp"

#include <cstdio>
#include <cstring>
#include <exception>
#include <iostream>
#include <iomanip>

std::ofstream Logging::_file{};

const std::array<std::string,3> Logging::TYPE_LABELS = {"INFO L", "WARN L", "ERROR L"};

std::mutex Logging::_lock{};

void Logging::log(LOG_TYPE type, int level, const std::string& message)
{
    std::lock_guard<std::mutex> lock(_lock);
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    _file << TYPE_LABELS[type] << level << " [" << std::put_time(&tm, "%a %b %d %H:%M:%S %Y") << "]: " << message << std::endl;
}

void Logging::debug_log(int level, const std::string& message)
{
    #ifdef DEBUG
        std::lock_guard<std::mutex> lock(_lock);
        std::time_t t = std::time(nullptr);
        std::tm tm = *std::localtime(&t);
        _file << "DEBUG L" << level << " [" << std::put_time(&tm, "%a %b %d %H:%M:%S %Y") << "]: " << message << std::endl;
    #endif
}


namespace output
{

void write_track_csv(const std::string& filename, const Track& track)
{
    std::ofstream file(filename);
    file << std::setprecision(5);
    for(auto i = 0u; i < track.t.size(); i++)
    {
        file << track.t[i];
        for(auto x : track.x[i])
        {
            file << " " << x;
        }
        file <<"\n";
    }
    file << std::flush;
}

template<>
void write_recording_csv<14>(const std::string& filename, const Recording<14>& recording)
{
    std::ofstream file(filename);
    file << std::setprecision(10);
    for(auto i = 0u; i < recording.times.size(); i++)
    {
        file << recording.times[i];
        for(auto j = 0u; j < 14; j++)
        {
            file << " " << recording.states[i][j];
        }
        file <<"\n";
    }
    file << std::flush;
}

void write_simulation_result(const std::string& filename, const Simulation::SimulationResult& result)
{
    FILE* file = fopen(filename.c_str(),"w");
    if(!file)
    {
        std::cout << "couldn't open file at: " << filename << std::endl;
    }
    const char* FORMAT = "%8.1f %+12.3f %+12.3f %+12.3f %+12.6f %+12.6f %+12.6f %+12.10f %+12.10f %+12.10f %+12.10f %+8.5f %+8.5f %+8.5f %11.3f %8.5f %+9.7f %+8.5f %+10.2f";
    for(auto i = 0u; i < result.times.size(); i++)
    {
        const auto& _x = result.ecef_states[i].x;
        const auto& _d = result.derived[i];
        
        fprintf(file, FORMAT, result.times[i], _x[0], _x[1], _x[2], _x[3], _x[4], _x[5], _x[6], _x[7], _x[8], _x[9],
             _x[10], _x[11], _x[12], _x[13], _d.mach, _d.flight_angle, _d.g_force, _d.altitude );

        for(auto val : result.other_data[i])
        {
            fprintf(file, " %+12.6f", val);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

} // end namespace