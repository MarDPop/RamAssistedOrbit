
#include "constants.hpp"
#include "output.hpp"
#include "simulation.hpp"
#include "test.hpp"

#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    std::string config_file = "./cfg/config.json";
    if(argc > 1)
    {
        if(strcmp(argv[1], "-t") == 0)
        {
            Tests::runAll();
            return 1;
        }
        config_file = std::string(argv[1]);
    }

    freopen("./output/run.log", "w", stdout);

    Logging::init();

    auto config = Simulation::load_config(config_file);
    auto result = Simulation::run(config);
    std::cout << "Writing output..." << std::endl;
    output::write_simulation_result(config["OUTPUT_DIRECTORY"].template get<std::string>() + "result.dat", result);
}