#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"

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
    std::vector<uint8_t> in_collision;              // whether or not particles are in collision - use uint8 instead of bool to avoid parallel writes to the same byte 

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit ParticlePool(unsigned capacity)
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
        , in_collision(capacity)
    {

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

        // assume particle is not initially in collision
        in_collision[slot] = 0;

        return slot;
    }

    /** Adds a new particle to the pool, with default initialized values for its properties. 
     * @returns the index of the particle in the pool
    */
    unsigned addParticle()
    {
        return addParticle(Vec3r::Zero(), 0);
    }

    /** Removes a particle from the pool.
     * @param slot : the index of the particle to remove
     */
    void removeParticle(unsigned slot)
    {
        freeSlot(slot);
    } 

};