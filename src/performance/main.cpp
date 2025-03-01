#include <iostream>
#include <fstream>
#include <string.h>
#include <string>

void show_help() {

}

struct run_options {
    int type;
    int fuel;
    double altitude = 0.0;
    double mach = 0.0;
    double throat_area = 1.0;
    double exit_area = 1.0;
    double inlet_area = 1.0;
    double nozzle_efficiency = 1.0;
    double combustion_efficiency = 1.0;
    double combustion_temperature = 2000.0;
};

run_options load_file(std::string filename) {

}

int main(int argc, char **argv) {
    if(argc < 1) {
        show_help();
        return -1;
    }

    if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        show_help();
        return 0;
    }

    if(strcmp(argv[1], "-f") == 0) {
        
    }
}