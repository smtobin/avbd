#pragma once

#include "common/common.hpp"
#include "simulation/SimulationContext.hpp"

#include "common/cuda/DeviceParticlePool.cuh"
#include "common/cuda/DeviceStaticParticleAdjacency.cuh"
#include "common/cuda/DeviceColorList.cuh"
#include "energy/cuda/DeviceEnergyRegistry.cuh"

namespace Sim
{


struct SimulationContextView
{
    ParticlePoolView particles;
    Energy::EnergyRegistryView energies;

    DeviceStaticParticleAdjacencyView static_adjacency;
    DeviceColorListView coloring;

    /** TODO: (10/28/26) remove this once we have a way to get device num colors to host  */
    unsigned num_colors;
};

////////////////////////////////////////////////////////////////////////////////////////////////

struct DeviceSimulationContext
{
    // memory pools
    DeviceParticlePool particles;
    Energy::DeviceEnergyRegistry energies;

    // adjacency information
    DeviceStaticParticleAdjacency static_adjacency;
    DeviceColorList coloring;

    /** Store number of colors
     * TODO: (10/28/26) find better way of getting device number of colors back to host
     * We can get away with keeping a static value for now because the scene is static
     */
    unsigned num_colors;

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

    SimulationContextView view()
    {
        SimulationContextView v;
        v.particles = particles.view();
        v.energies = energies.view();
        v.static_adjacency = static_adjacency.view();
        v.coloring = coloring.view();

        v.num_colors = num_colors;  // remove
        return v;
    }

    void upload(const SimulationContext& src)
    {
        particles.upload(src.particles);
        energies.upload(src.energies);
        static_adjacency.upload(src.static_adjacency);
        coloring.upload(src.coloring);

        num_colors = src.coloring.num_colors;  // remove

        params = src.params;
    }
};

} // namespace Sim