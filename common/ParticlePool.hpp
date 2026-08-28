#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"

/** Extra information required for oriented particles */
struct RotationPool : TombstonePool
{
    std::vector<Quaternion> rotations;              // particle rotations
    std::vector<Quaternion> buffered_rotations;     // particle buffered rotations
    std::vector<Quaternion> inertial_rotations;     // particle inertial rotations
    std::vector<Quaternion> previous_rotations;     // particle previous rotations
    std::vector<Quaternion> last_iter_rotations;    // particle positions at the end of the previous iteration (useful for Chebyshev acceleration)
    std::vector<Quaternion> last_last_iter_rotations;   // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    std::vector<Vec3r> angular_velocities;     // particle angular velocities
    std::vector<Vec3r> previous_angular_velocities;    // particle previous angular velocities
    std::vector<Vec3r> rotational_inertias;         // (body-frame, diagonal) particle rotational inertia

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit RotationPool(unsigned capacity)
        : TombstonePool(capacity)
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

    /** Adds new rotational information to the pool. */
    unsigned addParticle(const Quaternion& rotation, const Vec3r& rotational_inertia)
    {
        unsigned slot = allocSlot();
        rotations[slot] = rotation;
        inertial_rotations[slot] = rotation;
        previous_rotations[slot] = rotation;
        last_iter_rotations[slot] = rotation;
        last_last_iter_rotations[slot] = rotation;

        angular_velocities[slot] = Vec3r::Zero();
        previous_angular_velocities[slot] = Vec3r::Zero();
        rotational_inertias[slot] = rotational_inertia;

        return slot;
    }

    /** Adds a new particle to the pool, with default initialized values for its properties. 
     * @returns the index of the particle in the pool
    */
    unsigned addParticle()
    {
        return addParticle(Quaternion::Identity(), Vec3r::Zero());
    }

    /** Removes a particle from the pool.
     * @param slot : the index of the particle to remove
     */
    void removeParticle(unsigned slot)
    {
        freeSlot(slot);
    } 
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/** Defines a memory pool for all particles in the sim. */
struct ParticlePool : TombstonePool
{
    std::vector<Vec3r> positions;                   // particle positions
    std::vector<Vec3r> buffered_positions;          // position buffer to avoid race conditions for 
    std::vector<Vec3r> inertial_positions;          // particle inertial positions ('y' in the VBD paper)
    std::vector<Vec3r> previous_positions;          // particle previous positions
    std::vector<Vec3r> last_iter_positions;         // particle position at the end of the previous iteration (useful for Chebyshev acceleration)
    std::vector<Vec3r> last_last_iter_positions;    // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    std::vector<Vec3r> velocities;                  // particle velocities
    std::vector<Vec3r> previous_velocities;         // particle previous velocities
    std::vector<Real> masses;                       // particle masses
    std::vector<uint8_t> fixed;                     // whether or not the particle is fixed

    std::vector<unsigned> rotation_idx;             // index of the particle in the rotation pool. UINT_MAX if not oriented
    RotationPool rotation_pool;

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit ParticlePool(unsigned capacity, unsigned oriented_capacity)
        : TombstonePool(capacity)
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
        return rotation_pool.rotations[rotation_idx[p_idx]];
    }
    Quaternion& rotation(unsigned p_idx)
    {
        return rotation_pool.rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the buffered rotation of an oriented particle */
    const Quaternion& bufferedRotation(unsigned p_idx) const
    {
        return rotation_pool.buffered_rotations[rotation_idx[p_idx]];
    }
    Quaternion& bufferedRotation(unsigned p_idx)
    {
        return rotation_pool.buffered_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the previous rotation of an oriented particle */
    const Quaternion& previousRotation(unsigned p_idx) const
    {
        return rotation_pool.previous_rotations[rotation_idx[p_idx]];
    }
    Quaternion& previousRotation(unsigned p_idx)
    {
        return rotation_pool.previous_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get the angular velocity of an oriented particle */
    const Vec3r& angularVelocity(unsigned p_idx) const
    {
        return rotation_pool.angular_velocities[rotation_idx[p_idx]];
    }
    Vec3r& angularVelocity(unsigned p_idx)
    {
        return rotation_pool.angular_velocities[rotation_idx[p_idx]];
    }

    /** Convenience function to get the previous angular velocity of an oriented particle */
    const Vec3r& previousAngularVelocity(unsigned p_idx) const
    {
        return rotation_pool.previous_angular_velocities[rotation_idx[p_idx]];
    }
    Vec3r& previousAngularVelocity(unsigned p_idx)
    {
        return rotation_pool.previous_angular_velocities[rotation_idx[p_idx]];
    }

    /** Convenience function to get the rotational inertia of an oriented particle */
    const Vec3r& rotationalInertia(unsigned p_idx) const
    {
        return rotation_pool.rotational_inertias[rotation_idx[p_idx]];
    }
    Vec3r& rotationalInertia(unsigned p_idx)
    {
        return rotation_pool.rotational_inertias[rotation_idx[p_idx]];
    }

    /** Convenience function to get inertial rotation of an oriented particle */
    const Quaternion& inertialRotation(unsigned p_idx) const
    {
        return rotation_pool.inertial_rotations[rotation_idx[p_idx]];
    }
    Quaternion& inertialRotation(unsigned p_idx)
    {
        return rotation_pool.inertial_rotations[rotation_idx[p_idx]];
    }

    /** Convenience function to get previous iteration rotations of an oriented particle */
    const Quaternion& lastIterRotation(unsigned p_idx) const
    {
        return rotation_pool.last_iter_rotations[rotation_idx[p_idx]];
    }
    Quaternion& lastIterRotation(unsigned p_idx)
    {
        return rotation_pool.last_iter_rotations[rotation_idx[p_idx]];
    }
    const Quaternion& lastLastIterRotation(unsigned p_idx) const
    {
        return rotation_pool.last_last_iter_rotations[rotation_idx[p_idx]];
    }
    Quaternion& lastLastIterRotation(unsigned p_idx)
    {
        return rotation_pool.last_last_iter_rotations[rotation_idx[p_idx]];
    }

    /** Adds a new particle to the pool, given an initial position and mass.
     * Velocities are initialized to 0, and the particle is assumed to not be in collision.
     * 
     * @param position : the initial position of the particle
     * @param mass : the mass of the particle
     * @returns the index of the particle in the pool
     */
    unsigned addParticle(const Vec3r& position, Real mass)
    {
        unsigned slot = allocSlot();
        positions[slot] = position;
        masses[slot] = mass;

        // initialize other positional quantities to be the initial position
        inertial_positions[slot] = position;
        previous_positions[slot] = position;
        last_iter_positions[slot] = position;
        last_last_iter_positions[slot] = position;

        // initialize velocities
        velocities[slot] = Vec3r::Zero();
        previous_velocities[slot] = Vec3r::Zero();

        // assume particle is not fixed
        fixed[slot] = false;

        // particle is not oriented, so its rotation idx is UINT_MAX
        rotation_idx[slot] = std::numeric_limits<unsigned>::max();

        return slot;
    }

    /** Adds a new particle to the pool, with default initialized values for its properties. 
     * @returns the index of the particle in the pool
    */
    unsigned addParticle()
    {
        return addParticle(Vec3r::Zero(), 0);
    }

    /** Adds a new oriented particle to the pool, given an initial position, rotation, mass, and rotational inertia.
     * Velocities and angular velocities are initialized to 0, and particle is assumed to not be in collision.
     * @param position initial position of new particle
     * @param rotation initial rotation of new particle
     * @param mass translational mass of new particle
     * @param rotational_inertia diagonal, body-frame particle rotational inertia, expressed as a 3-vector
     */
    unsigned addOrientedParticle(const Vec3r& position, const Quaternion& rotation, Real mass, const Vec3r& rotational_inertia)
    {
        unsigned slot = addParticle(position, mass);
        rotation_idx[slot] = rotation_pool.addParticle(rotation, rotational_inertia);
        
        return slot;
    }

    /** Adds a new oriented particle to the pool, with default initialized values for its properties.
     * @returns the index of the particle in the pool
     */
    unsigned addOrientedParticle()
    {
        return addOrientedParticle(Vec3r::Zero(), Quaternion::Identity(), 0, Vec3r::Zero());
    }

    /** Removes a particle from the pool.
     * @param slot : the index of the particle to remove
     */
    void removeParticle(unsigned slot)
    {
        // if particle is oriented, free its slot in the rotational pool
        if (isOriented(slot))
            rotation_pool.removeParticle(rotation_idx[slot]);

        freeSlot(slot);
    } 

};

