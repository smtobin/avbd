#pragma once

#include "common/common.hpp"
#include "energy/EnergyRegistry.hpp"
#include "common/ParticlePool.hpp"

#include <vector>
#include <array>
#include <unordered_set>

template <typename T>
using PerEnergy = std::array<T, (unsigned)EnergyType::size>;

struct ParticleAdjacency {
    // For each particle, the range [adj_offsets[v], adj_offsets[v+1])
    // gives the list of constraint references in adj_entries
    std::vector<PerEnergy<unsigned>> e_offsets;   // size: numParticles + 1
    std::vector<PerEnergy<unsigned>> e_valences;

    struct Entry {
        EnergyType energy_type;  // type of energy
        unsigned energy_idx;   // index into that type's pool
        unsigned short  local_vertex_idx;  // the index of this vertex in the energy
    };
    std::vector<Entry> e_entries;  // size: sum of valences across all vertices

    /** Particle-particle adjacency */
    std::vector<unsigned> p_offsets;    // CSR offsets
    std::vector<unsigned> p_neighbors;  // indices of neighbors
    std::vector<unsigned> p_descending_valence;  // lists particle indices by decreasing valence

    /** Build vertex-energy adjacency and vertex-vertex adjacency. */
    void buildAdjacency(const ParticlePool& particle_pool, const Energy::EnergyRegistry& energy_registry)
    {
        unsigned num_particles = particle_pool.highest_index + 1;
        e_valences.resize(num_particles);
        e_valences.assign(num_particles, {});

        // step 1: count valence of each particle
        std::vector<std::unordered_set<unsigned>> adj_p_set(num_particles);

        energy_registry.forEachEnergyType([&] (const auto& pool) {
            using Pool = base_type_t<decltype(pool)>;
            unsigned e_type_ind = (unsigned)Pool::Type;

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

                        adj_p_set[p_idx].insert(pool.data[e_idx].particle_indices[k2]);
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
            for (unsigned e = 0; e < (unsigned)EnergyType::size; e++)
            {
                if ( e == (unsigned)EnergyType::size - 1 )
                {
                    e_offsets[i+1][0] = e_offsets[i][e] + e_valences[i][e];
                }
                else
                {
                    e_offsets[i][e+1] = e_offsets[i][e] + e_valences[i][e];
                }
            }
            
            p_offsets[i+1] = p_offsets[i] + adj_p_set[i].size();
        }

        // step 3: scatter entries for each energy
        e_entries.resize(e_offsets.back()[0]);  // the first entry in the last offsets array is the final size
        p_neighbors.resize(p_offsets.back());
        std::vector<PerEnergy<unsigned>> cursor = e_offsets;     // cursor per particle for where to put the next entry
        energy_registry.forEachEnergyType([&] (const auto& pool) {
            using Pool = base_type_t<decltype(pool)>;
            unsigned e_type_ind = (unsigned)Pool::Type;

            // iterate through each (active) energy in the pool
            for (unsigned e_idx : pool)
            {
                // iterate through each particle in the energy
                for (unsigned short k = 0; k < pool.NumParticlesPerEnergy; k++)
                {
                    unsigned p_idx = pool.data[e_idx].particle_indices[k];
                    e_entries[cursor[p_idx][e_type_ind]++] = {pool.Type, e_idx, k};
                }
            }
        });

        // step 4: vertex adjacency
        for (unsigned p_idx : particle_pool)
        {
            unsigned w = p_offsets[p_idx];
            for (const auto& v : adj_p_set[p_idx])
            {
                p_neighbors[w++] = v;
            }
        }

        // step 5: order particles by valence
        p_descending_valence.resize(num_particles);
        orderParticlesByDescendingValence();
    }

    /** Sort particles by valence. Useful for better coloring. */
    void orderParticlesByDescendingValence()
    {
        std::iota(p_descending_valence.begin(), p_descending_valence.end(), 0);
        std::sort(p_descending_valence.begin(), p_descending_valence.end(), [&](unsigned a, unsigned b) {
            unsigned deg_a = p_offsets[a+1] - p_offsets[a];
            unsigned deg_b = p_offsets[b+1] - p_offsets[b];
            return deg_a > deg_b;
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