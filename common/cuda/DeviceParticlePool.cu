#include "common/ParticlePool.hpp"
#include "common/cuda/DeviceParticlePool.cuh"

void DeviceParticlePool::upload(const ParticlePool& src)
{
    if (src.capacity != capacity)
        throw std::runtime_error("Host and device pool capacities do not match");

    positions.upload(src.positions, capacity);
    buffered_positions.upload(src.buffered_positions, capacity);
    inertial_positions.upload(src.inertial_positions, capacity);
    previous_positions.upload(src.previous_positions, capacity);
    last_iter_positions.upload(src.last_iter_positions, capacity);
    last_last_iter_positions.upload(src.last_last_iter_positions, capacity);
    velocities.upload(src.velocities, capacity);
    previous_velocities.upload(src.previous_velocities, capacity);
    masses.upload(src.masses, capacity);
    fixed.upload(src.fixed, capacity);
    rotation_idx.upload(src.rotation_idx, capacity);

    rotation_pool.upload(src.rotation_pool);
}