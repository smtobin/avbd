#pragma once

#include "common/common.hpp"
#include "common/ParticleAdjacency.hpp"

#include <vector>

struct ColorList
{
    static constexpr unsigned UNCOLORED = std::numeric_limits<unsigned>::max();

    unsigned num_colors;                // number of colors
    std::vector<unsigned> color;       // per-particle color
    std::vector<unsigned> color_counts;     // number of particles per color
    std::vector<unsigned> color_offsets;    // per-color offsets
    std::vector<unsigned> work_list;    // particle indices, grouped by color

    /** Greedy coloring of the particles given the adjacency structure.
     * Also sets the number of colors in the coloring.
     */
    void greedyColor(const ParticleAdjacency& adjacency, unsigned num_particles)
    {
        // resize the colors to the number of particles
        color.resize(num_particles);
        // reset all colors
        color.assign(num_particles, UNCOLORED);

        // scratch space for storing which colors are "used" by adjacent vertices
        std::vector<bool> used_color;
        for (unsigned v : adjacency.p_descending_valence)
        {
            std::cout << " Coloring vertex " << v << std::endl;
            std::cout << "    valence: " << adjacency.p_offsets[v+1] - adjacency.p_offsets[v] << std::endl;
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

    void buildColorList(const ParticleAdjacency& adjacency, unsigned num_particles)
    {
        // first, perform greedy coloring on the particle adjacency structure
        greedyColor(adjacency, num_particles);

        // count each color
        color_counts.resize(num_colors);
        color_counts.assign(num_colors, 0);
        for (unsigned c : color)
        {
            std::cout << "Color: " << c << std::endl;
            color_counts[c]++;
        }

        // use counts to generate offsets
        color_offsets.resize(num_colors+1, 0);
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
};