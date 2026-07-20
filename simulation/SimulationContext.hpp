#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/OrientedParticlePool.hpp"
#include "common/ParticleAdjacency.hpp"
#include "common/ColorList.hpp"
#include "energy/EnergyRegistry.hpp"
#include "collision/CollisionPrimitivePool.hpp"
#include "collision/LBVH.hpp"
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

    // collision detection
    Collision::CollisionPrimitivePool collision_pool;
    Collision::LBVH lbvh;

    // simulation parameters
    SimulationParams params;

    SimulationContext()
        : particles(1000)
        , oriented_particles(1000)
        , energies(1000)
        , collision_pool(1000, 1000)
    {
        
    }

    SimulationContext(
        unsigned particles_capacity,
        unsigned oriented_particles_capacity, 
        unsigned energies_capacity, 
        unsigned collision_primitive_capacity,
        unsigned collision_sdf_capacity
    )
     : particles(particles_capacity)
     , oriented_particles(oriented_particles_capacity)
     , energies(energies_capacity)
     , collision_pool(collision_primitive_capacity, collision_sdf_capacity)
    {
        
    }
};

} // namespace Sim