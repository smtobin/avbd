#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "common/ParticlePool.hpp"

/** Defines a memory pool for all particles in the sim. */
struct OrientedParticlePool : ParticlePool
{
    std::vector<Quaternion> rotations;              // particle rotations
    std::vector<Quaternion> inertial_rotations;     // particle inertial rotations
    std::vector<Quaternion> previous_rotations;     // particle previous rotations
    std::vector<Quaternion> last_iter_rotations;    // particle positions at the end of the previous iteration (useful for Chebyshev acceleration)
    std::vector<Quaternion> last_last_iter_rotations;   // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    std::vector<Vec3r> angular_velocities;     // particle angular velocities
    std::vector<Vec3r> previous_angular_velocities;    // particle previous angular velocities
    std::vector<Vec3r> rotational_inertias;         // (body-frame, diagonal) particle rotational inertia
    std::vector<uint8_t> fixed;                     // whether or not the particle is fixed

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit OrientedParticlePool(unsigned capacity)
        : ParticlePool(capacity)
        , rotations(capacity)
        , inertial_rotations(capacity)
        , previous_rotations(capacity)
        , last_iter_rotations(capacity)
        , last_last_iter_rotations(capacity)
        , angular_velocities(capacity)
        , previous_angular_velocities(capacity)
        , rotational_inertias(capacity)
        , fixed(capacity)
    {

    }

    /** Adds a new particle to the pool, given an initial position and mass.
     * Velocities are initialized to 0, and the particle is assumed to not be in collision.
     * 
     * @param position : the initial position of the particle
     * @param mass : the mass of the particle
     * @returns the index of the particle in the pool
     */
    unsigned addParticle(const Vec3r& position, const Quaternion& rot, Real mass, const Vec3r& I)
    {
        unsigned slot = ParticlePool::addParticle(position, mass);
        rotations[slot] = rot;
        rotational_inertias[slot] = I;

        inertial_rotations[slot] = rot;
        previous_rotations[slot] = rot;
        last_iter_rotations[slot] = rot;
        last_last_iter_rotations[slot] = rot;
        angular_velocities[slot] = Vec3r::Zero();
        previous_angular_velocities[slot] = Vec3r::Zero();
        fixed[slot] = 0;

        return slot;
    }

    /** Adds a new particle to the pool, with default initialized values for its properties. 
     * @returns the index of the particle in the pool
    */
    unsigned addParticle()
    {
        return addParticle(Vec3r::Zero(), Quaternion::Identity(), 0, Vec3r::Zero());
    }

    /** Removes a particle from the pool.
     * @param slot : the index of the particle to remove
     */
    void removeParticle(unsigned slot)
    {
        ParticlePool::removeParticle(slot);
    } 

};