#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "common/ParticlePool.hpp"
#include "common/ParticlePoolView.hpp"

// forward declarations of CPU types
struct RotationPool;
struct ParticlePool;

struct DeviceRotationPool : DeviceTombstonePool
{
    DeviceBuffer<Quaternion> rotations;
    DeviceBuffer<Quaternion> buffered_rotations;
    DeviceBuffer<Quaternion> inertial_rotations;
    DeviceBuffer<Quaternion> previous_rotations;
    DeviceBuffer<Quaternion> last_iter_rotations;
    DeviceBuffer<Quaternion> last_last_iter_rotations;
    DeviceBuffer<Vec3r> angular_velocities;
    DeviceBuffer<Vec3r> previous_angular_velocities;
    DeviceBuffer<Vec3r> rotational_inertias;

    explicit DeviceRotationPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , rotations(capacity)
        , buffered_rotations(capacity)
        , inertial_rotations(capacity)
        , previous_rotations(capacity)
        , last_iter_rotations(capacity)
        , last_last_iter_rotations(capacity)
        , angular_velocities(capacity)
        , previous_angular_velocities(capacity)
        , rotational_inertias(capacity)
    {
    }

    /** Fill out device view */
    RotationPoolView view()
    {
        RotationPoolView v;

        // tombstone pool data
        DeviceTombstonePool::view(v);

        // rotation pool data
        v.rotations = rotations.data();
        v.buffered_rotations = buffered_rotations.data();
        v.inertial_rotations = inertial_rotations.data();
        v.previous_rotations = previous_rotations.data();
        v.last_iter_rotations = last_iter_rotations.data();
        v.last_last_iter_rotations = last_last_iter_rotations.data();
        v.angular_velocities = angular_velocities.data();
        v.previous_angular_velocities = previous_angular_velocities.data();
        v.rotational_inertias = rotational_inertias.data();

        return v;
    }

    /** Upload data from host to device */
    void upload(const RotationPool& src)
    {
        if (src.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");

        DeviceTombstonePool::upload(src);

        rotations.upload(src.rotations.data(), capacity);
        buffered_rotations.upload(src.buffered_rotations.data(), capacity);
        inertial_rotations.upload(src.inertial_rotations.data(), capacity);
        previous_rotations.upload(src.previous_rotations.data(), capacity);
        last_iter_rotations.upload(src.last_iter_rotations.data(), capacity);
        last_last_iter_rotations.upload(src.last_last_iter_rotations.data(), capacity);
        angular_velocities.upload(src.angular_velocities.data(), capacity);
        previous_angular_velocities.upload(src.previous_angular_velocities.data(), capacity);
        rotational_inertias.upload(src.rotational_inertias.data(), capacity);
    }

    /** Synchronous download from device to host */
    void download(RotationPool& dst) const
    {
        if (dst.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");
        
        DeviceTombstonePool::download(dst);

        rotations.download(dst.rotations.data(), capacity);
        buffered_rotations.download(dst.buffered_rotations.data(), capacity);
        inertial_rotations.download(dst.inertial_rotations.data(), capacity);
        previous_rotations.download(dst.previous_rotations.data(), capacity);
        last_iter_rotations.download(dst.last_iter_rotations.data(), capacity);
        last_last_iter_rotations.download(dst.last_last_iter_rotations.data(), capacity);
        angular_velocities.download(dst.angular_velocities.data(), capacity);
        previous_angular_velocities.download(dst.previous_angular_velocities.data(), capacity);
        rotational_inertias.download(dst.rotational_inertias.data(), capacity);
    }
};

//////////////////////////////////////////////////////////////////////

struct DeviceParticlePool : DeviceTombstonePool
{
    DeviceBuffer<Vec3r> positions;                   // particle positions
    DeviceBuffer<Vec3r> buffered_positions;          // position buffer to avoid race conditions for 
    DeviceBuffer<Vec3r> inertial_positions;          // particle inertial positions ('y' in the VBD paper)
    DeviceBuffer<Vec3r> previous_positions;          // particle previous positions
    DeviceBuffer<Vec3r> last_iter_positions;         // particle position at the end of the previous iteration (useful for Chebyshev acceleration)
    DeviceBuffer<Vec3r> last_last_iter_positions;    // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    DeviceBuffer<Vec3r> velocities;                  // particle velocities
    DeviceBuffer<Vec3r> previous_velocities;         // particle previous velocities
    DeviceBuffer<Real> masses;                       // particle masses
    DeviceBuffer<uint8_t> fixed;                     // whether or not the particle is fixed

    DeviceBuffer<unsigned> rotation_idx;             // index of the particle in the rotation pool. UINT_MAX if not oriented
    DeviceRotationPool rotation_pool;

    explicit DeviceParticlePool(unsigned capacity, unsigned oriented_capacity)
        : DeviceTombstonePool(capacity)
        , positions(capacity)
        , buffered_positions(capacity)
        , inertial_positions(capacity)
        , previous_positions(capacity)
        , last_iter_positions(capacity)
        , last_last_iter_positions(capacity)
        , velocities(capacity)
        , previous_velocities(capacity)
        , masses(capacity)
        , fixed(capacity)
        , rotation_idx(capacity)
        , rotation_pool(oriented_capacity)
    {

    }

    /** Fill out device view */
    ParticlePoolView view()
    {
        ParticlePoolView v;

        // tombstone pool data
        DeviceTombstonePool::view(v);

        // particle pool data
        v.positions = positions.data();
        v.buffered_positions = buffered_positions.data();
        v.inertial_positions = inertial_positions.data();
        v.previous_positions = previous_positions.data();
        v.last_iter_positions = last_iter_positions.data();
        v.last_last_iter_positions = last_last_iter_positions.data();
        v.velocities = velocities.data();
        v.previous_velocities = previous_velocities.data();
        v.masses = masses.data();
        v.fixed = fixed.data();
        
        v.rotation_idx = rotation_idx.data();

        // rotation pool data
        v.rotation_pool_view = rotation_pool.view();

        return v;
    }

    /** Upload data from host to device */
    void upload(const ParticlePool& src)
    {
        if (src.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");

        DeviceTombstonePool::upload(src);

        positions.upload(src.positions.data(), capacity);
        buffered_positions.upload(src.buffered_positions.data(), capacity);
        inertial_positions.upload(src.inertial_positions.data(), capacity);
        previous_positions.upload(src.previous_positions.data(), capacity);
        last_iter_positions.upload(src.last_iter_positions.data(), capacity);
        last_last_iter_positions.upload(src.last_last_iter_positions.data(), capacity);
        velocities.upload(src.velocities.data(), capacity);
        previous_velocities.upload(src.previous_velocities.data(), capacity);
        masses.upload(src.masses.data(), capacity);
        fixed.upload(src.fixed.data(), capacity);
        rotation_idx.upload(src.rotation_idx.data(), capacity);

        rotation_pool.upload(src.rotation_pool);
    }

    /** Synchronous download from device to host */
    void download(ParticlePool& dst) const
    {
        if (dst.capacity != capacity)
            throw std::runtime_error("Host and device pool capacities do not match");

        DeviceTombstonePool::download(dst);

        positions.download(dst.positions.data(), capacity);
        buffered_positions.download(dst.buffered_positions.data(), capacity);
        inertial_positions.download(dst.inertial_positions.data(), capacity);
        previous_positions.download(dst.previous_positions.data(), capacity);
        last_iter_positions.download(dst.last_iter_positions.data(), capacity);
        last_last_iter_positions.download(dst.last_last_iter_positions.data(), capacity);
        velocities.download(dst.velocities.data(), capacity);
        previous_velocities.download(dst.previous_velocities.data(), capacity);
        masses.download(dst.masses.data(), capacity);
        fixed.download(dst.fixed.data(), capacity);
        rotation_idx.download(dst.rotation_idx.data(), capacity);
        
        rotation_pool.download(dst.rotation_pool);
    }
};