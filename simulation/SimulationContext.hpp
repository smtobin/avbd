#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/ParticleAdjacency.hpp"
#include "energy/EnergyRegistry.hpp"

namespace Simulation
{

/** Storage of all state in the sim. */
struct SimulationContext
{
    // memory pools
    ParticlePool particles;
    Energy::EnergyRegistry energies;

    // adjacency information for particles
    ParticleAdjacency adjacency;

    SimulationContext(unsigned capacity)
     : particles(capacity)
     , energies(capacity)
    {
        
    }
};

} // namespace Simulation