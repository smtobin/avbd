#pragma once

#include "common/common.hpp"

#include <vector>

struct ColorList
{
    std::vector<unsigned> color_offsets;    // per-color offsets
    std::vector<unsigned> work_list;    // particle indices, grouped by color

    /** TODO: greedy coloring given particle adjacency */
};