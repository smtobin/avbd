#pragma once

#include "common/ColorList.hpp"
#include "common/cuda/DeviceBuffer.cuh"
#include "common/cuda/DeviceValue.cuh"


struct DeviceColorListView
{
    unsigned* num_colors;

    unsigned* color;
    unsigned* color_counts;
    unsigned* color_offsets;
    unsigned* work_list;
};

//////////////////////////////////////////////////////////////////

struct DeviceColorList
{
    static constexpr unsigned UNCOLORED = ColorList::UNCOLORED;

    DeviceValue<unsigned> num_colors;       // number of active colors
    DeviceBuffer<unsigned> color;           // per-particle color
    DeviceBuffer<unsigned> color_counts;    // numer of particles per color
    DeviceBuffer<unsigned> color_offsets;   // per-color offsets
    DeviceBuffer<unsigned> work_list;       // particle indices, grouped by color

    /** TODO: (10/19/26) Port over conflicts, incremental recoloring. */

    DeviceColorList(unsigned particles_capacity)
        : num_colors(0u)
        , color(particles_capacity)
    {

    }

    DeviceColorListView view() const
    {
        DeviceColorListView v;
        v.num_colors = num_colors.data();
        v.color = color.data();
        v.color_counts = color_counts.data();
        v.color_offsets = color_offsets.data();
        v.work_list = color.data();
    }

    void upload(const ColorList& src)
    {
        // allocate memory
        color_counts.allocate(src.color_counts.size());
        color_offsets.allocate(src.color_offsets.size());
        work_list.allocate(src.work_list.size());

        // copy memory
        num_colors.set(src.num_colors);
        color.upload(src.color.data(), src.color.size());
        color_counts.upload(src.color_counts.data(), src.color_counts.size());
        color_offsets.upload(src.color_offsets.data(), src.color_offsets.size());
        work_list.upload(src.work_list.data(), src.work_list.size());
    }
};