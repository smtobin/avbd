#pragma once

#include "common/TombstonePool.hpp"

namespace Energy
{

/** TODO: switch to AoS */
struct HardConstraintEnergyInfo
{
    Real k;         // the finite stiffnesses of the quadratic energy
    Real lambda;    // the Lagrange multipliers enforcing the constraints
    Real C_prev;    // the constraint violation at the end of the previous time step
};

/** Pool of memory for HardConstraintEnergies.
 * Base class for all energies associated with hard constraints.
 * Really just has k (finite stiffness) and lambda (lagrange multipliers)
 */
template <typename EnergyInfo>
struct HardConstraintEnergyPool : TombstonePool
{
    Real k_start;           // the initial stiffness
    Real lambda_min;        // lower lambda bound for this constraint (if applicable)
    Real lambda_max;         // upper lambda bound for this constraint (if applicable)
    std::vector<EnergyInfo> data;

    explicit HardConstraintEnergyPool(
        unsigned capacity,
        Real k_start_,
        Real lambda_min_=std::numeric_limits<Real>::lowest(),
        Real lambda_max_=std::numeric_limits<Real>::max()
    )
        : TombstonePool(capacity)
        , k_start(k_start_)
        , lambda_min(lambda_min_)
        , lambda_max(lambda_max_)
        , data(capacity)
    {}

    /** Add an energy
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy()
    {
        unsigned slot = allocSlot();

        // initialize lambda and k
        data[slot].lambda = 0;
        data[slot].k = k_start;
        data[slot].C_prev = 0;

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