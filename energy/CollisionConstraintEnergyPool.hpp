#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"

namespace Energy
{

struct CollisionConstraintEnergyInfo
{
    Real k_n; Real k_t; Real k_b;        // the finite stiffnesses of the quadratic energy (one for each direction of the collision basis)
    Real lambda_n; Real lambda_t; Real lambda_b;   // the Lagrange multipliers enforcing the constraints (normal, tangent, binormal)
    Real C_n_prev; Real C_t_prev; Real C_b_prev;   // the constraint violation at the end of the previous time step

    /** The collision orthonormal frame */
    Vec3r normal;   // collision normal vector
    Vec3r tangent;  // collision tangent vector
    Vec3r binormal; // collision binormal vector
    
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
    unsigned addEnergy()
    {
        unsigned slot = allocSlot();

        // initialize lambda and k
        data[slot].lambda_n = 0;
        data[slot].lambda_t = 0;
        data[slot].lambda_b = 0;
        data[slot].k_n = 0;
        data[slot].k_t = 0;
        data[slot].k_b = 0;
        data[slot].C_n_prev = 0;
        data[slot].C_t_prev = 0;
        data[slot].C_b_prev = 0;

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