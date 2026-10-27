#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "energy/collision/CollisionConstraintEnergyDefs.hpp"

namespace Energy
{

/** Pool of memory specifically for collision-related constraint energies */
template <typename CollisionEnergyInfo>
struct CollisionConstraintEnergyPool : TombstonePool
{
    std::vector<CollisionEnergyInfo> data;

    explicit CollisionConstraintEnergyPool(
        unsigned capacity
    )
        : TombstonePool(capacity)
        , data(capacity)
    {}

    /** Returns a non-owning view into the collision constraint energy pool */
    void view(CollisionConstraintEnergyPoolView<CollisionEnergyInfo>& v) const
    {
        TombstonePool::view(v);
        v.data() = data.data();
    }

    /** Add an energy
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(Real k_start, Real mu_s, Real mu_k)
    {
        unsigned slot = allocSlot();

        // initialize lambdas and ks
        data[slot].lambda_n = 0;
        data[slot].lambda_t = 0;
        data[slot].lambda_b = 0;
        data[slot].k_start = k_start;
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