#include "simulation/Simulation.hpp"

#include <gmsh.h>


int main(int argc, char **argv) 
{
    gmsh::initialize();

    if (argc > 1)
    {
        std::string config_filename(argv[1]);
        Config::SimulationConfig config(YAML::LoadFile(config_filename));
        Sim::Simulation sim(config);
        return sim.run();
    }
    else
    {
        std::cerr << "No config file specified!" << std::endl;
    }
}