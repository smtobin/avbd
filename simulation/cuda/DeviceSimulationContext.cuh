#pragma once

#include "common/common.hpp"
#include "simulation/SimulationContext.hpp"

#include "common/cuda/DeviceParticlePool.cuh"
#include "common/cuda/DeviceStaticParticleAdjacency.cuh"
#include "common/cuda/DeviceColorList.cuh"
#include "energy/cuda/DeviceEnergyRegistry.cuh"

namespace Sim
{

struct DeviceSimulationContext
{
    // memory pools
    DeviceParticlePool particles;
    Energy::DeviceEnergyRegistry energies;

    // adjacency information
    DeviceStaticParticleAdjacency static_adjacency;
    DeviceColorList coloring;

    /** TODO: (10/23/26) Collision detection things */

    // simulation parameters
    SimulationParams params;

    DeviceSimulationContext(
        unsigned particles_capacity,
        unsigned oriented_particles_capacity,
        unsigned energies_capacity,
        unsigned collision_primitive_capacity,
        unsigned collision_sdf_capacity
    )
        : particles(particles_capacity, oriented_particles_capacity)
        , energies(energies_capacity)
        , coloring(particles_capacity)
    {

    }

    void upload(const SimulationContext& src)
    {
        particles.upload(src.particles);
        energies.upload(src.energies);
        static_adjacency.upload(src.static_adjacency);
        coloring.upload(src.coloring);

        params = src.params;
    }
};

} // namespace Sim