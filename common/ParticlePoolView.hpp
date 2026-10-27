#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"

struct RotationPoolView : TombstonePoolView
{
    Quaternion* rotations;
    Quaternion* buffered_rotations;
    Quaternion* inertial_rotations;
    Quaternion* previous_rotations;
    Quaternion* last_iter_rotations;
    Quaternion* last_last_iter_rotations;
    Vec3r* angular_velocities;
    Vec3r* previous_angular_velocities;
    Vec3r* rotational_inertias;
};

struct ParticlePoolView : TombstonePoolView
{
    Vec3r* positions;
    Vec3r* buffered_positions;
    Vec3r* inertial_positions;
    Vec3r* previous_positions;
    Vec3r* last_iter_positions;
    Vec3r* last_last_iter_positions;
    Vec3r* velocities;
    Vec3r* previous_velocities;
    float* masses;
    uint8_t* fixed;

    unsigned* rotation_idx;
    RotationPoolView rotation_pool_view;

    /** Queries if particle is oriented or not
     * @param p_idx the index of the particle in the pool
     * @returns if the particle is oriented
     */
    bool isOriented(unsigned p_idx) const 
    {
        return rotation_idx[p_idx] != std::numeric_limits<unsigned>::max();
    }

    /** Convenience function to get the rotation of an oriented particle */
    const Quaternion& rotation(unsigned p_idx) const
    {
        return rotation_pool_view.rotations[rotation_idx[p_idx]];
    }
    Quaternion& rotation(unsigned p_idx)
    {
        return rotation_pool_view.rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the buffered rotation of an oriented particle */
    const Quaternion& bufferedRotation(unsigned p_idx) const
    {
        return rotation_pool_view.buffered_rotations[rotation_idx[p_idx]];
    }
    Quaternion& bufferedRotation(unsigned p_idx)
    {
        return rotation_pool_view.buffered_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the previous rotation of an oriented particle */
    const Quaternion& previousRotation(unsigned p_idx) const
    {
        return rotation_pool_view.previous_rotations[rotation_idx[p_idx]];
    }
    Quaternion& previousRotation(unsigned p_idx)
    {
        return rotation_pool_view.previous_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the angular velocity of an oriented particle */
    const Vec3r& angularVelocity(unsigned p_idx) const
    {
        return rotation_pool_view.angular_velocities[rotation_idx[p_idx]];
    }
    Vec3r& angularVelocity(unsigned p_idx)
    {
        return rotation_pool_view.angular_velocities[rotation_idx[p_idx]];
    }

    /** Convenience function to get the previous angular velocity of an oriented particle */
    const Vec3r& previousAngularVelocity(unsigned p_idx) const
    {
        return rotation_pool_view.previous_angular_velocities[rotation_idx[p_idx]];
    }
    Vec3r& previousAngularVelocity(unsigned p_idx)
    {
        return rotation_pool_view.previous_angular_velocities[rotation_idx[p_idx]];
    }

    /** Convenience function to get the rotational inertia of an oriented particle */
    const Vec3r& rotationalInertia(unsigned p_idx) const
    {
        return rotation_pool_view.rotational_inertias[rotation_idx[p_idx]];
    }
    Vec3r& rotationalInertia(unsigned p_idx)
    {
        return rotation_pool_view.rotational_inertias[rotation_idx[p_idx]];
    }

    /** Convenience function to get inertial rotation of an oriented particle */
    const Quaternion& inertialRotation(unsigned p_idx) const
    {
        return rotation_pool_view.inertial_rotations[rotation_idx[p_idx]];
    }
    Quaternion& inertialRotation(unsigned p_idx)
    {
        return rotation_pool_view.inertial_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get previous iteration rotations of an oriented particle */
    const Quaternion& lastIterRotation(unsigned p_idx) const
    {
        return rotation_pool_view.last_iter_rotations[rotation_idx[p_idx]];
    }
    Quaternion& lastIterRotation(unsigned p_idx)
    {
        return rotation_pool_view.last_iter_rotations[rotation_idx[p_idx]];
    }
    const Quaternion& lastLastIterRotation(unsigned p_idx) const
    {
        return rotation_pool_view.last_last_iter_rotations[rotation_idx[p_idx]];
    }
    Quaternion& lastLastIterRotation(unsigned p_idx)
    {
        return rotation_pool_view.last_last_iter_rotations[rotation_idx[p_idx]];
    }
};