#pragma once
 
#include "common/common.hpp"
#include "energy/EnergyRegistry.hpp"
#include "common/ParticlePool.hpp"

#include <vector>
#include <array>

template <typename T>
using PerDynamicEnergy = std::array<T, (unsigned)DynamicEnergyType::count>;

/** Structure for tracking particle-energy and particle-particle adjacency from collision constraints.
 * Collision constraints change every frame, so keep this separate from the main particle adjacency, which only needs a rebuild if topology changes.
 */
struct DynamicParticleAdjacency
{
    std::vector<PerDynamicEnergy<unsigned>> e_offsets;
    std::vector<PerDynamicEnergy<unsigned>> e_valences;

    struct Entry {
        DynamicEnergyType energy_type;  // type of energy
        unsigned energy_idx;   // index into that type's pool
        unsigned short  local_vertex_idx;  // the index of this vertex in the energy
    };
    std::vector<Entry> e_entries;  // size: sum of valences across all vertices

    /** Particle-particle adjacency */
    std::vector<unsigned> p_offsets;    // CSR offsets
    std::vector<unsigned> p_degree;    // number of neighbors per particle (include duplicates)
    std::vector<unsigned> p_neighbors;  // indices of neighbors (may have duplicates!)

    void buildAdjacency(const ParticlePool& particle_pool, const Energy::EnergyRegistry& energy_registry)
    {
        unsigned num_particles = particle_pool.totalSize();
        e_valences.assign(num_particles, {});
        p_degree.assign(num_particles, 0);

        // step 1: count valence of each particle
        energy_registry.forEachDynamicEnergyType([&] (const auto& pool) {
            using Pool = base_type_t<decltype(pool)>;
            unsigned e_type_ind = (unsigned)Pool::DynamicType;

            // iterate through each (active) energy in the pool
            for (unsigned e_idx : pool)
            {
                // iterate through each particle in the energy
                for (unsigned k = 0; k < pool.NumParticlesPerEnergy; k++)
                {
                    // adjacent energies
                    unsigned p_idx = pool.data[e_idx].particle_indices[k];
                    e_valences[p_idx][e_type_ind]++;

                    // adjacent vertices
                    for (unsigned k2 = 0; k2 < pool.NumParticlesPerEnergy; k2++)
                    {
                        if (k == k2)
                            continue;

                        p_degree[p_idx]++;
                    }
                }
            }
        });

        // step 2: prefix sum to get offsets
        e_offsets.resize(num_particles + 1);
        e_offsets[0][0] = 0;
        p_offsets.resize(num_particles + 1);
        p_offsets[0] = 0;
        for (unsigned i = 0; i < num_particles; i++)
        {
            // std::cout << "valence " << i << ": " << e_valences[i] << std::endl;
            for (unsigned e = 0; e < (unsigned)DynamicEnergyType::count; e++)
            {
                if ( e == (unsigned)DynamicEnergyType::count - 1 )
                {
                    e_offsets[i+1][0] = e_offsets[i][e] + e_valences[i][e];
                }
                else
                {
                    e_offsets[i][e+1] = e_offsets[i][e] + e_valences[i][e];
                }
            }
            
            p_offsets[i+1] = p_offsets[i] + p_degree[i];
        }

        // step 3: scatter entries for each energy
        e_entries.resize(e_offsets.back()[0]);
        p_neighbors.resize(p_offsets.back());
        std::vector<PerDynamicEnergy<unsigned>> e_cursor = e_offsets;
        std::vector<unsigned> p_cursor = p_offsets;
        energy_registry.forEachDynamicEnergyType([&] (const auto& pool) {
            using Pool = base_type_t<decltype(pool)>;
            unsigned e_type_ind = (unsigned)Pool::DynamicType;

            // iterate through each (active) energy in the pool
            for (unsigned e_idx : pool)
            {
                // iterate through each particle in the energy
                for (unsigned short k = 0; k < pool.NumParticlesPerEnergy; k++)
                {
                    unsigned p_idx = pool.data[e_idx].particle_indices[k];
                    e_entries[e_cursor[p_idx][e_type_ind]++] = {pool.DynamicType, e_idx, k};

                    // adjacent vertices
                    for (unsigned k2 = 0; k2 < pool.NumParticlesPerEnergy; k2++)
                    {
                        if (k == k2)
                            continue;

                        p_neighbors[p_cursor[p_idx]++] = pool.data[e_idx].particle_indices[k2];
                    }
                }
            }
        });
    }

};