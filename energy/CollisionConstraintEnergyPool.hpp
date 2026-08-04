#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"

namespace Energy
{

struct CollisionConstraintEnergyInfo
{
    Vec3r k;        // the finite stiffnesses of the quadratic energy (one for each direction of the collision basis)
    Vec3r lambda;   // the Lagrange multipliers enforcing the constraints (tangent, binormal, normal)
    Vec3r C_prev;   // the constraint violation at the end of the previous time step
    Vec3r tangent;  // 
    Vec3r binormal;
    Vec3r normal;
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
        data[slot].lambda = Vec3r::Zero();
        data[slot].k = Vec3r::Constant(k_start);
        data[slot].C_prev = Vec3r::Zero();

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