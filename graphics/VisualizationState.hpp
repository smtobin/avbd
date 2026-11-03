#pragma once

#include "common/common.hpp"

namespace Graphics
{

/** Stores the relevant information for visualization of the current sim state. */
struct VisualizationState
{
    unsigned* counts = nullptr;   // counts[0] = number of particles, counts[1] = number of oriented particles

    Vec3r* particle_positions = nullptr;
    unsigned* particle_rotation_idx = nullptr;
    Quaternion* particle_rotations = nullptr;
};

#ifdef __CUDACC__

void allocHost(VisualizationState& s, unsigned max_particles, unsigned max_oriented_particles)
{
    CUDA_CHECK(cudaMallocHost(&s.particle_positions,        max_particles * sizeof(Vec3r)));
    CUDA_CHECK(cudaMallocHost(&s.particle_rotation_idx,     max_particles * sizeof(unsigned)));
    CUDA_CHECK(cudaMallocHost(&s.particle_rotations,        max_oriented_particles * sizeof(Quaternion)));
    CUDA_CHECK(cudaMallocHost(&s.counts,                    2 * sizeof(unsigned)));
}

void allocDevice(VisualizationState& s, size_t max_points, size_t max_oriented)
{
    CUDA_CHECK(cudaMalloc(&s.particle_positions,    max_points * sizeof(Vec3r)));
    CUDA_CHECK(cudaMalloc(&s.particle_rotation_idx, max_particles * sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&s.particle_rotations,    max_oriented_particles * sizeof(Quaternion)));
    CUDA_CHECK(cudaMalloc(&s.counts,                2 * sizeof(unsigned)));
}

#endif

} // namespace Graphics
