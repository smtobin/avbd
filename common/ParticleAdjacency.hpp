#pragma once

#include "common/common.hpp"

#include <vector>

struct ParticleAdjacency {
    // For each vertex, the range [adj_offsets[v], adj_offsets[v+1])
    // gives the list of constraint references in adj_entries
    std::vector<unsigned> adj_offsets;   // size: numParticles + 1
    
    struct Entry {
        EnergyType energy_type;  // type of energy
        unsigned energy_idx;   // index into that type's pool
        unsigned short  local_vertex_idx;  // the index of this vertex in the energy
    };
    std::vector<Entry> adj_entries;  // size: sum of valences across all vertices

    void addEntry(unsigned p_idx, EnergyType e_type, unsigned e_idx, unsigned short local_idx)
    {
        // get the current end of the range for this particle
        unsigned end = adj_offsets[p_idx+1];

        // insert an entry at this spot
        adj_entries.insert(end, {e_type, e_idx, local_idx});

        // increment the end of the range for this particle and all particles following
        for (unsigned idx = p_idx+1; idx < adj_offsets.size(); idx++)
            adj_offsets[idx]++;
    }

    template <typename EnergyPool>
    void buildAdjacency(const EnergyPool& pool)
    {
        for (unsigned e_idx = 0; e_idx < pool.highest_index; e_idx++)
        {
            // particle indices affected by the energy
            const auto& indices = pool.particle_indices[c_idx];

            for (unsigned k = 0; k < indices.size(); k++)
            {
                addEntry(indices[k], EnergyPool::Type, e_idx, k);
            }
        }
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