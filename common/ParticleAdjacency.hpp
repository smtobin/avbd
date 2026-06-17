#pragma once

#include "common/common.hpp"
#include "energy/EnergyRegistry.hpp"
#include "common/ParticlePool.hpp"

#include <vector>

struct ParticleAdjacency {
    // For each vertex, the range [adj_offsets[v], adj_offsets[v+1])
    // gives the list of constraint references in adj_entries
    std::vector<unsigned> adj_offsets;   // size: numParticles + 1
    std::vector<unsigned> valences;

    struct Entry {
        EnergyType energy_type;  // type of energy
        unsigned energy_idx;   // index into that type's pool
        unsigned short  local_vertex_idx;  // the index of this vertex in the energy
    };
    std::vector<Entry> adj_entries;  // size: sum of valences across all vertices

    void buildAdjacency(const ParticlePool& particle_pool, const Energy::EnergyRegistry& energy_registry)
    {
        unsigned num_particles = particle_pool.highest_index;
        valences.resize(num_particles);
        valences.assign(num_particles, 0);

        // step 1: count valence of each particle
        energy_registry.forEachEnergyType([&] (const auto& pool) {
            // iterate through each (active) energy in the pool
            for (unsigned e_idx : pool)
            {
                // iterate through each particle in the energy
                for (unsigned k = 0; k < pool.NumParticlesPerEnergy; k++)
                {
                    unsigned p_idx = pool.particle_indices[e_idx][k];
                    valences[p_idx]++;
                }
            }
        });

        // step 2: prefix sum to get offsets
        adj_offsets.resize(num_particles + 1);
        adj_offsets[0] = 0;
        for (unsigned i = 0; i < num_particles; i++)
        {
            adj_offsets[i+1] = adj_offsets[i] + valences[i];
        }

        // step 3: scatter entries for each energy
        std::vector<unsigned> cursor = adj_offsets;     // cursor per particle for where to put the next entry
        energy_registry.forEachEnergyType([&] (const auto& pool) {
            // iterate through each (active) energy in the pool
            for (unsigned e_idx : pool)
            {
                // iterate through each particle in the energy
                for (unsigned short k = 0; k < pool.NumParticlesPerEnergy; k++)
                {
                    unsigned p_idx = pool.particle_indices[e_idx][k];
                    adj_entries[cursor[p_idx]++] = {pool.Type, e_idx, k};
                }
            }
        });
    }
};

/**
uint32_t start = adj.adjOffsets[pIdx];
uint32_t end   = adj.adjOffsets[pIdx + 1];

for (uint32_t e = start; e < end; e++) {
    auto& entry = adj.adjEntries[e];
    if (entry.constraintType == STRETCH) {
        // read stretchConstraints[entry.constraintIdx]
        // entry.localVertexIdx tells you whether you're i or j
    } else if (entry.constraintType == VOLUME) {
        // ...
    }
}
 */