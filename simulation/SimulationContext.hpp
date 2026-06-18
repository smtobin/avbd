#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/ParticleAdjacency.hpp"
#include "energy/EnergyRegistry.hpp"

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

    SimulationContext()
        : particles(1000)
        , energies(1000)
    {
        
    }

    SimulationContext(unsigned particles_capacity, unsigned energies_capacity)
     : particles(particles_capacity)
     , energies(energies_capacity)
    {
        
    }
};

} // namespace Sim