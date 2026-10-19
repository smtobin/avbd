#pragma once

#include "common/cuda/DeviceParticlePool.cuh"
#include "common/StaticParticleAdjacency.hpp"

struct DeviceStaticParticleAdjacencyView
{
    using Entry = StaticParticleAdjacency::Entry;

    PerStaticEnergy<unsigned>* e_offsets;
    PerStaticEnergy<unsigned>* e_valences;
    Entry* e_entries;

    unsigned* p_offsets;
    unsigned* p_degree;
    unsigned* p_neighbors;
    unsigned* p_descending_valence;
};

////////////////////////////////////////////////////////////////////

struct DeviceStaticParticleAdjacency
{
    using Entry = StaticParticleAdjacency::Entry;

    DeviceBuffer<PerStaticEnergy<unsigned>> e_offsets;
    DeviceBuffer<PerStaticEnergy<unsigned>> e_valences;
    DeviceBuffer<Entry> e_entries;

    DeviceBuffer<unsigned> p_offsets;
    DeviceBuffer<unsigned> p_degree;
    DeviceBuffer<unsigned> p_neighbors;
    DeviceBuffer<unsigned> p_descending_valence;

    DeviceStaticParticleAdjacencyView view() const
    {
        DeviceStaticParticleAdjacencyView v;

        v.e_offsets = e_offsets.data();
        v.e_valences = e_valences.data();
        v.e_entries = e_entries.data();

        v.p_offsets = p_offsets.data();
        v.p_degree = p_degree.data();
        v.p_neighbors = p_neighbors.data();
        v.p_descending_valence = p_descending_valence.data();

        return v;
    }

    void upload(const StaticParticleAdjacency& src)
    {
        // allocate memory
        e_offsets.allocate(src.e_offsets.size());
        e_valences.allocate(src.e_valences.size());
        e_entries.allocate(src.e_entries.size());

        p_offsets.allocate(src.p_offsets.size());
        p_degree.allocate(src.p_degree.size());
        p_neighbors.allocate(src.p_neighbors.size());
        p_descending_valence.allocate(src.p_descending_valence.size());

        // copy data
        e_offsets.upload(src.e_offsets.data(), src.e_offsets.size());
        e_valences.upload(src.e_valences.data(), src.e_valences.size());
        e_entries.upload(src.e_entries.data(), src.e_entries.size());

        p_offsets.upload(src.p_offsets.data(), src.p_offsets.size());
        p_degree.upload(src.p_degree.data(), src.p_degree.size());
        p_neighbors.upload(src.p_neighbors.data(), src.p_neighbors.size());
        p_descending_valence.upload(src.p_descending_valence.data(), src.p_descending_valence.size());
    }

};