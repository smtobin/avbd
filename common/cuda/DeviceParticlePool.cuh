#pragma once

#include "common/cuda/Vec3.cuh"
#include "common/cuda/Quaternion.cuh"
#include "common/cuda/DeviceTombstonePool.cuh"
#include "common/ParticlePool.hpp"


struct DeviceRotationPoolView : DeviceTombstonePoolView
{
    DQuaternion* rotations;
    DQuaternion* buffered_rotations;
    DQuaternion* inertial_rotations;
    DQuaternion* previous_rotations;
    DQuaternion* last_iter_rotations;
    DQuaternion* last_last_iter_rotations;
    DVec3* angular_velocities;
    DVec3* previous_angular_velocities;
    DVec3* rotational_inertias;
};

struct DeviceParticlePoolView : DeviceTombstonePoolView
{
    DVec3* positions;
    DVec3* buffered_positions;
    DVec3* inertial_positions;
    DVec3* previous_positions;
    DVec3* last_iter_positions;
    DVec3* last_last_iter_positions;
    DVec3* velocities;
    DVec3* previous_velocities;
    float* masses;
    uint8_t* fixed;

    unsigned* rotation_idx;
    DeviceRotationPoolView rotation_pool_view;

    /** Queries if particle is oriented or not
     * @param p_idx the index of the particle in the pool
     * @returns if the particle is oriented
     */
    bool isOriented(unsigned p_idx) const 
    {
        return rotation_idx[p_idx] != std::numeric_limits<unsigned>::max();
    }

    /** Convenience function to get the rotation of an oriented particle */
    const DQuaternion& rotation(unsigned p_idx) const
    {
        return rotation_pool_view.rotations[rotation_idx[p_idx]];
    }
    DQuaternion& rotation(unsigned p_idx)
    {
        return rotation_pool_view.rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the buffered rotation of an oriented particle */
    const DQuaternion& bufferedRotation(unsigned p_idx) const
    {
        return rotation_pool_view.buffered_rotations[rotation_idx[p_idx]];
    }
    DQuaternion& bufferedRotation(unsigned p_idx)
    {
        return rotation_pool_view.buffered_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the previous rotation of an oriented particle */
    const DQuaternion& previousRotation(unsigned p_idx) const
    {
        return rotation_pool_view.previous_rotations[rotation_idx[p_idx]];
    }
    DQuaternion& previousRotation(unsigned p_idx)
    {
        return rotation_pool_view.previous_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the angular velocity of an oriented particle */
    const DVec3r& angularVelocity(unsigned p_idx) const
    {
        return rotation_pool_view.angular_velocities[rotation_idx[p_idx]];
    }
    DVec3r& angularVelocity(unsigned p_idx)
    {
        return rotation_pool_view.angular_velocities[rotation_idx[p_idx]];
    }

    /** Convenience function to get the previous angular velocity of an oriented particle */
    const DVec3r& previousAngularVelocity(unsigned p_idx) const
    {
        return rotation_pool_view.previous_angular_velocities[rotation_idx[p_idx]];
    }
    DVec3r& previousAngularVelocity(unsigned p_idx)
    {
        return rotation_pool_view.previous_angular_velocities[rotation_idx[p_idx]];
    }

    /** Convenience function to get the rotational inertia of an oriented particle */
    const DVec3r& rotationalInertia(unsigned p_idx) const
    {
        return rotation_pool_view.rotational_inertias[rotation_idx[p_idx]];
    }
    DVec3r& rotationalInertia(unsigned p_idx)
    {
        return rotation_pool_view.rotational_inertias[rotation_idx[p_idx]];
    }

    /** Convenience function to get inertial rotation of an oriented particle */
    const DQuaternion& inertialRotation(unsigned p_idx) const
    {
        return rotation_pool_view.inertial_rotations[rotation_idx[p_idx]];
    }
    DQuaternion& inertialRotation(unsigned p_idx)
    {
        return rotation_pool_view.inertial_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get previous iteration rotations of an oriented particle */
    const DQuaternion& lastIterRotation(unsigned p_idx) const
    {
        return rotation_pool_view.last_iter_rotations[rotation_idx[p_idx]];
    }
    DQuaternion& lastIterRotation(unsigned p_idx)
    {
        return rotation_pool_view.last_iter_rotations[rotation_idx[p_idx]];
    }
    const DQuaternion& lastLastIterRotation(unsigned p_idx) const
    {
        return rotation_pool_view.last_last_iter_rotations[rotation_idx[p_idx]];
    }
    DQuaternion& lastLastIterRotation(unsigned p_idx)
    {
        return rotation_pool_view.last_last_iter_rotations[rotation_idx[p_idx]];
    }
};

//////////////////////////////////////////////////////////////////////

struct DeviceRotationPool : DeviceTombstonePool
{
    DeviceBuffer<DQuaternion> rotations;
    DeviceBuffer<DQuaternion> buffered_rotations;
    DeviceBuffer<DQuaternion> inertial_rotations;
    DeviceBuffer<DQuaternion> previous_rotations;
    DeviceBuffer<DQuaternion> last_iter_rotations;
    DeviceBuffer<DQuaternion> last_last_iter_rotations;
    DeviceBuffer<DVec3> angular_velocities;
    DeviceBuffer<DVec3> previous_angular_velocities;
    DeviceBuffer<DVec3> rotational_inertias;

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
    DeviceRotationPoolView view() const
    {
        DeviceRotationPoolView v;

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

        rotations.upload(src.rotations, capacity);
        buffered_rotations.upload(src.buffered_rotations, capacity);
        inertial_rotations.upload(src.inertial_rotations, capacity);
        previous_rotations.upload(src.previous_rotations, capacity);
        last_iter_rotations.upload(src.last_iter_rotations, capacity);
        last_last_iter_rotations.upload(src.last_last_iter_rotations, capacity);
        angular_velocities.upload(src.angular_velocities, capacity);
        previous_angular_velocities.upload(src.previous_angular_velocities, capacity);
        rotational_inertias.upload(src.rotational_inertias, capacity);
    }
};

//////////////////////////////////////////////////////////////////////

struct DeviceParticlePool : DeviceTombstonePool
{
    std::vector<DVec3> positions;                   // particle positions
    std::vector<DVec3> buffered_positions;          // position buffer to avoid race conditions for 
    std::vector<DVec3> inertial_positions;          // particle inertial positions ('y' in the VBD paper)
    std::vector<DVec3> previous_positions;          // particle previous positions
    std::vector<DVec3> last_iter_positions;         // particle position at the end of the previous iteration (useful for Chebyshev acceleration)
    std::vector<DVec3> last_last_iter_positions;    // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    std::vector<DVec3> velocities;                  // particle velocities
    std::vector<DVec3> previous_velocities;         // particle previous velocities
    std::vector<float> masses;                       // particle masses
    std::vector<uint8_t> fixed;                     // whether or not the particle is fixed

    std::vector<unsigned> rotation_idx;             // index of the particle in the rotation pool. UINT_MAX if not oriented
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
    DeviceParticlePoolView view() const
    {
        DeviceParticlePoolView v;

        // tombstone pool data
        DeviceParticlePool::view(v);

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
};