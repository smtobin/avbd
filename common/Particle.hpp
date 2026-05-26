#pragma once

#include "common/common.hpp"
#include "common/Math.hpp"

/** Represents a general "particle" that has 3 DOF in a vector space (e.g. a particle with only positional DOF).
 * This particle can be a node in a volumteric continuum, or a purely positional node in a rod,
 *  or even some other expression of generalized coordinates.
 * 
 */
struct Particle
{
    constexpr static int DOF = 3;
    Vec3r position;         // position (global frame) of the particle
    Vec3r velocity;     // translational velocity (global frame) of the particle
    Real mass;              // mass of the particle
    Vec3r prev_position;    // previous position (at the end of the last time step) of the particle
    bool fixed=false;             // if true, the particle is "fixed" and should not move

    Particle(const Vec3r& position_, const Vec3r& velocity_, Real mass_, bool fixed_)
        : position(position_), velocity(velocity_), mass(mass_), prev_position(position_), fixed(fixed_)
    {
    }

    Particle(const Vec3r& position_)
        : position(position_), velocity(0,0,0), mass(1), prev_position(position_), fixed(false)
    {
    }

    /** Updates the particle based on its current velocity (in the absence of constraints) and applied external wrench.
     * @param dt - the time step
     * @param F_ext - external force
     */
    void inertialUpdate(Real dt, const Vec3r& F_ext);

    /** Updates the particle based on its current velocity (in the absence of constraints) with no applied external wrench.
     * @param dt - the time step
     */
    void inertialUpdate(Real dt);

    void inertialUpdateAcceleration(Real dt, const Vec3r& a);

    /** Updates the particle given some change in position and orientiation.
     * @param dpos - the change in position (specified in global coordinates)
     */
    void positionUpdate(const Vec3r& dpos);

    /** Updates the particle given some change in position and orientation.
     * @param dp - a 6-vector with the combined change in position (first 3 coords) and orientation (last 3 coords).
     */
    void positionUpdate(const Vec6r& dp);

    /** Updates the particle velocity given a time step, and previous position and orientation.
     * @param dt - the time step
     */
    void velocityUpdate(Real dt);
};