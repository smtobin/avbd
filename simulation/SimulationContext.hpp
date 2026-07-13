#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/OrientedParticlePool.hpp"
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
    OrientedParticlePool oriented_particles;
    Energy::EnergyRegistry energies;

    // adjacency information for particles
    ParticleAdjacency adjacency;
    ColorList coloring;

    // simulation parameters
    SimulationParams params;

    SimulationContext()
        : particles(1000)
        , oriented_particles(1000)
        , energies(1000)
    {
        
    }

    SimulationContext(unsigned particles_capacity, unsigned oriented_particles_capacity, unsigned energies_capacity)
     : particles(particles_capacity)
     , oriented_particles(oriented_particles_capacity)
     , energies(energies_capacity)
    {
        
    }
};

} // namespace Sim