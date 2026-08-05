#pragma once

#include "common/common.hpp"
#include "common/StaticParticleAdjacency.hpp"

#include <vector>

template <unsigned N>
struct ColorBitset
{
    uint64_t bitset[N];     // 64*N colors supported

    ColorBitset()
        : bitset({})
    {
    }

    inline void addColor(unsigned c)
    {
        bitset[c / 64] |= (1ull << (c % 64));
    }

    inline unsigned findFirstUnset()
    {
        for (unsigned w = 0; w < N; w++)
        {
            uint64_t inverted = ~bitset[w];
            if (inverted != 0)
                return w * 64  + std::countr_zero(inverted);
        }

        // every color in range is taken - throw an error
        throw std::runtime_error("Ran out of colors!");
    }
};

struct ColorList
{
    static constexpr unsigned UNCOLORED = std::numeric_limits<unsigned>::max();

    unsigned num_colors;                // number of colors
    std::vector<unsigned> color;       // per-particle color
    std::vector<unsigned> color_counts;     // number of particles per color
    std::vector<unsigned> color_offsets;    // per-color offsets
    std::vector<unsigned> work_list;    // particle indices, grouped by color

    std::vector<uint8_t> is_conflicted;         // tracks which particles have intra-color conflicts
    std::vector<unsigned> touched_conflicted;   // tracks which particles were "touched" by the conflict algorithm

    /** CSR structure for per-color conflicted particles */
    std::vector<unsigned> conflicted_by_color_offsets;  // size num_colors + 1
    std::vector<unsigned> conflicted_by_color_counts;   // per color conflicted counts, size = num_colors + 1
    std::vector<unsigned> conflicted_by_color_entries;  // size = leftover conflict count

    /** For incremental recoloring */
    std::vector<unsigned> dirty;    // list of particle indices that must be checked for conflicts. Populated by collision detection
    std::vector<unsigned> next_dirty;   // buffer of dirty vertices for next iteration of recoloring 
    std::vector<unsigned> last_dirty_round;     // the "round" of recoloring that this particle was last dirty
    std::vector<unsigned> candidate_color;      // buffered candidate colors for dirty vertices
    std::vector<unsigned> touched_this_frame;   // list of particles that were touched by the recoloring algorithm this frame

    static constexpr unsigned SENTINEL = 0xFFFFFFFFu;   // sentinel value for the last dirty round

    /** Initializes memory based on the max particle capacity in the particle pool */
    ColorList(unsigned particles_capacity)
        : color(particles_capacity, UNCOLORED)
        , is_conflicted(particles_capacity, 0)
        , last_dirty_round(particles_capacity, SENTINEL)
        , candidate_color(particles_capacity)
    {

    }

    /** Greedy coloring of the particles given the adjacency structure.
     * Also sets the number of colors in the coloring.
     * 
     * Uses Welsh-Powell greedy coloring - colors vertices in order according to their degree (number of adjacent vertices)
     */
    void greedyColor(const StaticParticleAdjacency& adjacency, unsigned num_particles)
    {
        // reset all colors within particle range
        std::fill(color.begin(), color.begin() + num_particles, UNCOLORED);

        // scratch space for storing which colors are "used" by adjacent vertices
        std::vector<bool> used_color;
        for (unsigned v : adjacency.p_descending_valence)
        {
            // std::cout << " Coloring vertex " << v << std::endl;
            // std::cout << "    valence: " << adjacency.p_offsets[v+1] - adjacency.p_offsets[v] << std::endl;
            // reset the used color buffer
            used_color.assign(used_color.size(), false);

            // iterate through adjacent vertices
            // mark colors as used in the scratch array
            for (unsigned e = adjacency.p_offsets[v]; e < adjacency.p_offsets[v+1]; e++)
            {
                unsigned neighbor = adjacency.p_neighbors[e];
                if (color[neighbor] != UNCOLORED)
                {
                    if (color[neighbor] >= used_color.size())
                    {
                        used_color.resize(color[neighbor]+1, false);
                    }
                    used_color[color[neighbor]] = true;
                }
            }

            // pick the lowest unused color
            unsigned c = 0;
            while (c < used_color.size() && used_color[c])
                c++;
            color[v] = c;
        }

        num_colors = used_color.size();
    }

    /** DSATUR greedy coloring
     * First sort particles by "saturation" (number of different adjacent colors),
     * with a tie-breaker of particle degree.
     * 
     * Supposedly produces a lower number of colors than Welsh-Powell.
     */
    void greedyColor2(const StaticParticleAdjacency& adjacency, unsigned num_particles)
    {
        std::fill(color.begin(), color.begin() + num_particles, UNCOLORED);

        std::vector<unsigned> saturation(num_particles, 0);

        // scratch space for storing which colors are "used" by adjacent vertices
        std::vector<bool> used_color;
        for (unsigned i = 0; i < num_particles; i++)
        {
            // find next particle to color
            unsigned best = UNCOLORED;
            for (unsigned j = 0; j < num_particles; j++)
            {
                // skip colored particles
                if (color[j] != UNCOLORED)
                    continue;
                
                // initially set the next particle to the first uncolored particle
                if (best == UNCOLORED)
                    best = j;
                
                // get the degree
                unsigned deg_best = adjacency.p_offsets[best+1] - adjacency.p_offsets[best];
                unsigned deg_j = adjacency.p_offsets[j+1] - adjacency.p_offsets[j];

                // check saturation
                if (saturation[j] > saturation[best])
                    best = j;
                
                // if tied on saturation, check degree
                else if (saturation[j] == saturation[best] && deg_j > deg_best)
                    best = j;
            }

            if (best == UNCOLORED)
                throw std::runtime_error("best color is UNCOLORED!");
                

            // reset the used color buffer
            used_color.assign(used_color.size(), false);

            // iterate through adjacent vertices
            // mark colors as used in the scratch array
            for (unsigned e = adjacency.p_offsets[best]; e < adjacency.p_offsets[best+1]; e++)
            {
                unsigned neighbor = adjacency.p_neighbors[e];
                if (color[neighbor] != UNCOLORED)
                {
                    if (color[neighbor] >= used_color.size())
                    {
                        used_color.resize(color[neighbor]+1, false);
                    }
                    used_color[color[neighbor]] = true;
                }
            }

            // pick the lowest unused color
            unsigned c = 0;
            while (c < used_color.size() && used_color[c])
                c++;
            color[best] = c;
            used_color.resize(std::max<unsigned>(used_color.size(), c+1));

            // update saturation of neighbors
            for (unsigned e = adjacency.p_offsets[best]; e < adjacency.p_offsets[best+1]; e++)
            {
                unsigned neighbor = adjacency.p_neighbors[e];
                if (color[neighbor] != UNCOLORED)
                    continue;

                used_color.assign(used_color.size(), false);
                for (unsigned k = adjacency.p_offsets[neighbor]; k < adjacency.p_offsets[neighbor+1]; k++)
                {
                    unsigned neighbor_neighbor = adjacency.p_neighbors[k];
                    if (color[neighbor_neighbor] != UNCOLORED)
                        used_color[color[neighbor_neighbor]] = true;
                }
                unsigned cnt = 0;
                for (auto used : used_color)
                {
                    if (used)
                        cnt++;
                }

                saturation[neighbor] = cnt;
            }
            
        }

        num_colors = used_color.size();

    }

    void countColors(unsigned num_particles)
    {
        // count each color
        color_counts.resize(num_colors);
        color_counts.assign(num_colors, 0);
        for (unsigned c_idx = 0; c_idx < num_particles; c_idx++)
        {
            // std::cout << "Color: " << c << std::endl;
            color_counts[color[c_idx]]++;
        }
    }

    /** Slight post-processing of colors
     * Lump small colors (<20% of largest color into other colors)
     */
    void mergeSmallColors(const StaticParticleAdjacency& adjacency, unsigned num_particles)
    {
        // find small colors (<20% of largest color)
        unsigned largest_color_size = 0;
        for (unsigned i = 0; i < num_colors; i++)
            if (color_counts[i] > largest_color_size)
                largest_color_size = color_counts[i];

        std::vector<unsigned> small_colors;
        for (unsigned i = 0; i < num_colors; i++)
        {
            if (color_counts[i] < largest_color_size/5)
                small_colors.push_back(i);
        }

        // simple strategy: for vertices in colors that are being merged, check adjacent colors and find the one with the least conflicts
        for (unsigned small_color : small_colors)
        {
            for (unsigned p_idx = 0; p_idx < num_particles; p_idx++)
            {
                unsigned c = color[p_idx];
                if (c == small_color)
                {
                    // look at neighbors and find color with fewest conflicts
                    unsigned best_color = UNCOLORED;
                    unsigned best_conflicts = std::numeric_limits<unsigned>::max();
                    for (unsigned candidate = 0; candidate < num_colors; candidate++)
                    {
                        unsigned conflicts = 0;
                        if (candidate == small_color)
                            continue;

                        for (unsigned e = adjacency.p_offsets[p_idx]; e < adjacency.p_offsets[p_idx+1]; e++)
                        {
                            unsigned neighbor = adjacency.p_neighbors[e];
                            if (color[neighbor] == candidate)
                                conflicts++;
                        }

                        if (conflicts < best_conflicts)
                        {
                            best_color = candidate;
                            best_conflicts = conflicts;
                        }
                    }
                    
                    // assign the color with least conflicts
                    color[p_idx] = best_color;
                    // mark this particle as dirty
                    dirty.push_back(p_idx);
                    last_dirty_round[p_idx] = 0;
                    touched_this_frame.push_back(p_idx);
                    
                }
            }
        }

        num_colors -= small_colors.size();
    }

    /** Builds an initial color list based on the static adjacency */
    void buildInitialColorList(const StaticParticleAdjacency& static_adjacency, const DynamicParticleAdjacency& dynamic_adjacency, unsigned num_particles)
    {
        // first, perform greedy coloring on the particle adjacency structure
        greedyColor2(static_adjacency, num_particles);

        /** TODO: (08/02/26) merge small colors? The commented code does not work properly. Idk if this something that actually is needed, however */
        // count each color (necessary for merging the small colors)
        // countColors(num_particles);

        // std::cout << "Initial coloring - " << num_colors << " colors: " << std::endl;
        // for (unsigned i = 0; i < num_colors; i++)
        //     std::cout << "  Color " << i << ": " << color_counts[i] << std::endl;

        
        // mergeSmallColors(static_adjacency, num_particles);

        // incrementalRecoloring(static_adjacency, dynamic_adjacency, num_particles);

        // std::cout << "After small color merging - " << num_colors << " colors: " << std::endl;
        // for (unsigned i = 0; i < num_colors; i++)
        //     std::cout << "  Color " << i << ": " << color_counts[i] << std::endl;

        rebuildWorkList(num_particles);
        
    }

    void rebuildWorkList(unsigned num_particles)
    {
        // recount colors
        countColors(num_particles);

        // use counts to generate offsets
        color_offsets.resize(num_colors+1, 0);
        color_offsets[0] = 0;
        for (unsigned c = 0; c < num_colors; c++)
            color_offsets[c+1] = color_offsets[c] + color_counts[c];

        // use offsets to generate work list
        work_list.resize(num_particles);
        std::vector<unsigned> cursor = color_offsets;
        for (unsigned v = 0; v < num_particles; v++)
        {
            work_list[cursor[color[v]]++] = v;
        }
    }

    /** Tries to mark a particle as dirty in the color list. If it has already been marked dirty this frame, it won't be added to the dirty list. */
    inline void markParticleDirty(unsigned p_idx)
    {
        // using last dirty round prevents duplication
        // i.e. dirty particles only added once
        // note: this assumes that last_dirty_round has enough allocated space for all particles in the sim
        if (last_dirty_round[p_idx] != 0)
        {
            last_dirty_round[p_idx] = 0;
            touched_this_frame.push_back(p_idx);
            dirty.push_back(p_idx);
        }
    }


    /** Incremental coloring to fix potential conflicts given the new particle adjacency */
    void incrementalRecoloring(const StaticParticleAdjacency& static_adj, const DynamicParticleAdjacency& dynamic_adj, unsigned num_particles)
    {
        constexpr unsigned MAX_COLOR_ITERS = 6;
        constexpr unsigned NUM_WORDS_IN_BITSET = 4;     // 4*64 = 256 colors supported

        // reset conflicted particles
        // these will still be marked "dirty" from last frame and checked again
        for (unsigned p_idx : touched_conflicted)
        {
            is_conflicted[p_idx] = 0;
        }
        touched_conflicted.clear();

        unsigned iter = 0;
        while (!dirty.empty() && iter < MAX_COLOR_ITERS)
        {
            for (unsigned p_idx : dirty)
            {
                // gather colors used by static and dynamic neighbors
                ColorBitset<NUM_WORDS_IN_BITSET> bitset;
                unsigned static_p_start = static_adj.p_offsets[p_idx];
                unsigned static_p_end = static_adj.p_offsets[p_idx+1];
                unsigned dyn_p_start = dynamic_adj.p_offsets[p_idx];
                unsigned dyn_p_end = dynamic_adj.p_offsets[p_idx+1];

                bool must_recolor = false;
                for (unsigned n_idx = static_p_start; n_idx < static_p_end; n_idx++)
                {
                    unsigned neighbor = static_adj.p_neighbors[n_idx];
                    unsigned c_neighbor = color[neighbor];
                    bitset.addColor(c_neighbor);
                    if (c_neighbor == color[p_idx] && neighbor < p_idx)
                        must_recolor = true;
                }

                for (unsigned n_idx = dyn_p_start; n_idx < dyn_p_end; n_idx++)
                {
                    unsigned neighbor = dynamic_adj.p_neighbors[n_idx];
                    unsigned c_neighbor = color[neighbor];
                    bitset.addColor(c_neighbor);
                    if (c_neighbor == color[p_idx] && neighbor < p_idx)
                        must_recolor = true;
                }

                // if there is a conflict and it is determined that this particle must recolor,
                // get the first unused color and use it as the candidate color
                if (must_recolor)
                {
                    candidate_color[p_idx] = bitset.findFirstUnset();
                    if (candidate_color[p_idx] >= num_colors)
                        num_colors++;
                }
                else
                {
                    candidate_color[p_idx] = color[p_idx];
                }
            }

            // apply and propagate - a changed color might create a new conflict with a neighbor that wasn't dirty before
            next_dirty.clear();
            for (unsigned p_idx : dirty)
            {
                if (candidate_color[p_idx] != color[p_idx])
                {
                    color[p_idx] = candidate_color[p_idx];

                    unsigned static_p_start = static_adj.p_offsets[p_idx];
                    unsigned static_p_end = static_adj.p_offsets[p_idx+1];
                    unsigned dyn_p_start = dynamic_adj.p_offsets[p_idx];
                    unsigned dyn_p_end = dynamic_adj.p_offsets[p_idx+1];

                    for (unsigned n_idx = static_p_start; n_idx < static_p_end; n_idx++)
                    {
                        unsigned neighbor = static_adj.p_neighbors[n_idx];
                        if (last_dirty_round[neighbor] != iter + 1)
                        {
                            last_dirty_round[neighbor] = iter + 1;
                            touched_this_frame.push_back(neighbor);
                            next_dirty.push_back(neighbor);
                        }
                        
                    }

                    for (unsigned n_idx = dyn_p_start; n_idx < dyn_p_end; n_idx++)
                    {
                        unsigned neighbor = dynamic_adj.p_neighbors[n_idx];
                        if (last_dirty_round[neighbor] != iter + 1)
                        {
                            last_dirty_round[neighbor] = iter + 1;
                            touched_this_frame.push_back(neighbor);
                            next_dirty.push_back(neighbor);
                        }
                    }
                }
            }

            dirty.swap(next_dirty);
            iter++;
        }

        // std::cout << "Remaining 'dirty' particles: " << dirty.size() << std::endl;

        // if there are leftover dirty particles, finalize the list of conflicts by checking the dirties against their neighbors
        if (!dirty.empty())
            finalizeConflicts(static_adj, dynamic_adj, dirty);
        // otherwise, we must still set up the empty conflict CSR structure
        else
        {
            conflicted_by_color_counts.assign(num_colors+1, 0);
            conflicted_by_color_offsets.assign(num_colors+1, 0);
        }

        // rebuild CSR
        rebuildWorkList(num_particles);

        // reset touched and last dirty lists
        // important to do this now because CollisionDetector will add dirty particles before recoloring is performed
        for (unsigned p : touched_this_frame)
        {
            last_dirty_round[p] = SENTINEL;
        }
        touched_this_frame.clear();

        // residual (still-unresolved) particles stay flagged at the seed marker so
        // _markParticleDirty won't re-queue them as duplicates before the next call
        for (unsigned p : dirty)
        {
            last_dirty_round[p] = 0;
        }
    }

    void finalizeConflicts(const StaticParticleAdjacency& static_adj, const DynamicParticleAdjacency& dynamic_adj, const std::vector<unsigned>& residual_dirty)
    {
        for (unsigned p_idx : residual_dirty)
        {
            /** TODO: (08/02/26) Replace with "forEachNeighbor" or some nice way to iterate over adjacency structure */
            unsigned static_p_start = static_adj.p_offsets[p_idx];
            unsigned static_p_end = static_adj.p_offsets[p_idx+1];
            unsigned dyn_p_start = dynamic_adj.p_offsets[p_idx];
            unsigned dyn_p_end = dynamic_adj.p_offsets[p_idx+1];

            for (unsigned n_idx = static_p_start; n_idx < static_p_end; n_idx++)
            {
                unsigned neighbor = static_adj.p_neighbors[n_idx];
                if (color[neighbor] == color[p_idx])
                {
                    if (!is_conflicted[p_idx])
                    {
                        is_conflicted[p_idx] = 1;
                        touched_conflicted.push_back(p_idx);
                    }
                    if (!is_conflicted[neighbor])
                    {
                        is_conflicted[neighbor] = 1;
                        touched_conflicted.push_back(neighbor);
                    }
                }
                
            }

            for (unsigned n_idx = dyn_p_start; n_idx < dyn_p_end; n_idx++)
            {
                unsigned neighbor = dynamic_adj.p_neighbors[n_idx];
                if (color[neighbor] == color[p_idx])
                {
                    if (!is_conflicted[p_idx])
                    {
                        is_conflicted[p_idx] = 1;
                        touched_conflicted.push_back(p_idx);
                    }
                    if (!is_conflicted[neighbor])
                    {
                        is_conflicted[neighbor] = 1;
                        touched_conflicted.push_back(neighbor);
                    }
                }
            }
        }

        /** Build CSR structure per color for in-conflict particles */
        // count number of conflicted particles in each color
        conflicted_by_color_counts.resize(num_colors);
        conflicted_by_color_counts.assign(num_colors, 0);
        for (unsigned p_idx : touched_conflicted)
        {
            conflicted_by_color_counts[color[p_idx]]++;
        }
        // use counts to generate offsets
        conflicted_by_color_offsets.resize(num_colors+1, 0);
        conflicted_by_color_offsets[0] = 0;
        for (unsigned c = 0; c < num_colors; c++)
        {
            conflicted_by_color_offsets[c+1] = conflicted_by_color_offsets[c] + conflicted_by_color_counts[c];
        }

        // use offsets to generate entries
        conflicted_by_color_entries.resize(conflicted_by_color_offsets.back());
        auto cursor = conflicted_by_color_offsets;
        for (unsigned p_idx : touched_conflicted)
        {
            conflicted_by_color_entries[cursor[color[p_idx]]++] = p_idx;
        }

    }
};