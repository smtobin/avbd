#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/ParticleAdjacency.hpp"
#include "common/ColorList.hpp"
#include "common/ThreadPool.hpp"
#include "energy/EnergyRegistry.hpp"
#include "simulation/SimulationParams.hpp"

namespace Sim
{

/** Storage of all state in the sim. */
struct SimulationContext
{
    // memory pools
    ParticlePool particles;
    Energy::EnergyRegistry energies;

    // adjacency information for particles
    ParticleAdjacency adjacency;
    ColorList coloring;

    // simulation parameters
    SimulationParams params;

    // thread pool
    ThreadPool thread_pool;

    SimulationContext()
        : particles(1000)
        , energies(1000)
        , thread_pool(4)
    {
        
    }

    SimulationContext(unsigned particles_capacity, unsigned energies_capacity)
     : particles(particles_capacity)
     , energies(energies_capacity)
     , thread_pool(4)
    {
        
    }
};

} // namespace Sim