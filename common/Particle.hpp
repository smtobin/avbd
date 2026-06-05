#pragma once

#include "common/common.hpp"
#include "common/Math.hpp"

#include <vector>

/** Represents a general "particle" that has 3 DOF in a vector space (e.g. a particle with only positional DOF).
 * This particle can be a node in a volumteric continuum, or a purely positional node in a rod,
 *  or even some other expression of generalized coordinates.
 * 
 */
struct Particle
{
    constexpr static int DOF = 3;
    Vec3r position;         // position (global frame) of the particle
    Vec3r inertial_position; // the inertially predicted position ('y' in the VBD paper)
    Vec3r velocity;     // translational velocity (global frame) of the particle
    Real mass;              // mass of the particle
    Vec3r prev_position;    // previous position (at the end of the last time step) of the particle
    Vec3r prev_velocity;    // previous velocity (from the end of two time steps ago)
    bool fixed=false;             // if true, the particle is "fixed" and should not move

    Real kd=0;          // coefficient of damping for this particle

    Vec3r last_iter_position;       // position at the end of the previous iteration
    Vec3r last_last_iter_position;  // position at the end of 2 iterations ago

    // whether the particle is actively in collision
    // particles in collision should be Chebyshev-acclerated
    // (mutable because we set this from inside energy evaluations, which are const -> TODO: is this a code smell?)
    mutable bool in_collision = false;  

    

    // energy expressions affecting this particle
    // stored as (energy ptr, index) pairs where the index is the index of this particle in the energy expression
    std::vector<std::pair<const Energy::Energy_Base*, int>> energies;

    Particle(const Vec3r& position_, const Vec3r& velocity_, Real mass_, bool fixed_, Real kd_)
        : position(position_), velocity(velocity_), mass(mass_),
         prev_position(position_), prev_velocity(velocity_), fixed(fixed_), kd(kd_),
         last_iter_position(position_), last_last_iter_position(position_),
         in_collision(false)
    {
    }

    Particle(const Vec3r& position_)
        : position(position_), velocity(0,0,0), mass(1),
         prev_position(position_), prev_velocity(velocity), fixed(false), kd(0),
         last_iter_position(position_), last_last_iter_position(position_),
         in_collision(false)
    {
    }

    /** Updates the particle based on its current velocity (in the absence of constraints) and an external acceleration.
     * @param dt - the time step
     * @param a_ext - external acceleration
     */
    void inertialUpdate(Real dt, const Vec3r& a_ext);

    /** Solves the individual particle system (the vertex block), and updates the particle position. */
    void solveParticle(Real dt);

    /** Updates the particle velocity given a time step, and previous position and orientation.
     * @param dt - the time step
     */
    void velocityUpdate(Real dt);
};