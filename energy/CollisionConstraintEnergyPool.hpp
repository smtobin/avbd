#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"

namespace Energy
{

struct CollisionConstraintEnergyInfo
{
    Real k_n; Real k_t; Real k_b;        // the finite stiffnesses of the quadratic energy (one for each direction of the collision basis)
    Real lambda_n; Real lambda_t; Real lambda_b;   // the Lagrange multipliers enforcing the constraints (normal, tangent, binormal)
    Real C_n_prev;   // the constraint violation at the end of the previous time step (just for the normal direction)
    // (we do not use the constraint error smoothing for the tangent or binormal directions because the forces are naturally capped by the friction cone)

    /** The collision orthonormal frame */
    Vec3r normal;   // collision normal vector
    Vec3r tangent;  // collision tangent vector
    Vec3r binormal; // collision binormal vector

    /** Friction properties */
    bool use_static; // whether or not static friction is active
    // tracks whether the friction Lagrange multiplier was unclamped during the last iteration 
    // if, at the end of the time step, this is true and we were using kinetic friction, we switch to static
    bool unclamped_lambda_tb;   
    Real mu_s;  // static friction coefficient
    Real mu_k;  // kinetic friction coefficient
    
};

/** Pool of memory specifically for collision-related constraint energies */
template <typename CollisionEnergyInfo>
struct CollisionConstraintEnergyPool : TombstonePool
{
    Real k_start;
    std::vector<CollisionEnergyInfo> data;

    explicit CollisionConstraintEnergyPool(
        unsigned capacity,
        Real k_start_
    )
        : TombstonePool(capacity)
        , k_start(k_start_)
        , data(capacity)
    {}

    /** Add an energy
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(Real mu_s, Real mu_k)
    {
        unsigned slot = allocSlot();

        // initialize lambdas and ks
        data[slot].lambda_n = 0;
        data[slot].lambda_t = 0;
        data[slot].lambda_b = 0;
        data[slot].k_n = k_start;
        data[slot].k_t = k_start;
        data[slot].k_b = k_start;
        // initialize C_prev
        data[slot].C_n_prev = 0;
        // initialize friciton coeffs
        data[slot].mu_s = mu_s;
        data[slot].mu_k = mu_k;
        data[slot].use_static = false;  // default to kinetic friction initially
        data[slot].unclamped_lambda_tb = false;

        return slot;
    }

    /** Remove an energy
     * @param slot : the index of the energy in the pool to remove
     */
    void removeEnergy(unsigned slot)
    {
        freeSlot(slot);
    }
        
};

} // namespace Energy