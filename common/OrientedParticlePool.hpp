#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "common/Quaternion.hpp"

/** Defines a memory pool for all particles in the sim. */
struct OrientedParticlePool : TombstonePool
{
    std::vector<Vec3r> positions;                   // particle positions
    std::vector<Quaternion> rotations;              // particle rotations
    std::vector<Vec3r> inertial_positions;          // particle inertial positions ('y' in the VBD paper)
    std::vector<Quaternion> inertial_rotations;     // particle inertial rotations
    std::vector<Vec3r> previous_positions;          // particle previous positions
    std::vector<Quaternion> previous_rotations;     // particle previous rotations
    std::vector<Vec3r> last_iter_positions;         // particle position at the end of the previous iteration (useful for Chebyshev acceleration)
    std::vector<Vec3r> last_last_iter_positions;    // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    std::vector<Quaternion> last_iter_rotations;    // particle positions at the end of the previous iteration (useful for Chebyshev acceleration)
    std::vector<Quaternion> last_last_iter_rotations;   // particle positions at the end of 2 iterations ago (useful for Chebyshev acceleration)
    std::vector<Vec3r> velocities;                  // particle velocities
    std::vector<Quaternion> angular_velocities;     // particle angular velocities
    std::vector<Vec3r> previous_velocities;         // particle previous velocities
    std::vector<Quaternion> previous_angular_velocities;    // particle previous angular velocities
    std::vector<Real> masses;                       // particle masses
    std::vector<Vec3r> rotational_inertias;         // (body-frame, diagonal) particle rotational inertia
    std::vector<uint8_t> in_collision;              // whether or not particles are in collision - use uint8 instead of bool to avoid parallel writes to the same byte 
    std::vector<uint8_t> fixed;                     // whether or not the particle is fixed

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit OrientedParticlePool(unsigned capacity)
        : TombstonePool(capacity)
        , positions(capacity)
        , rotations(capacity)
        , inertial_positions(capacity)
        , inertial_rotations(capacity)
        , previous_positions(capacity)
        , previous_rotations(capacity)
        , last_iter_positions(capacity)
        , last_last_iter_positions(capacity)
        , last_iter_rotations(capacity)
        , last_last_iter_rotations(capacity)
        , velocities(capacity)
        , angular_velocities(capacity)
        , previous_velocities(capacity)
        , previous_angular_velocities(capacity)
        , masses(capacity)
        , rotational_inertias(capacity)
        , in_collision(capacity)
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
        unsigned slot = allocSlot();
        positions[slot] = position;
        rotations[slot] = rot;
        masses[slot] = mass;
        rotational_inertias[slot] = I;


        // initialize other positional quantities to be the initial position
        inertial_positions[slot] = position;
        previous_positions[slot] = position;
        last_iter_positions[slot] = position;
        last_last_iter_positions[slot] = position;

        // initialize velocities
        velocities[slot] = Vec3r::Zero();
        previous_velocities[slot] = Vec3r::Zero();

        // assume particle is not initially in collision
        in_collision[slot] = 0;

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
        freeSlot(slot);
    } 

};