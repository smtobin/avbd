#pragma once

#include "common/TombstonePool.hpp"

namespace Energy
{

/** Pool of memory for HardConstraintEnergies.
 * Base class for all energies associated with hard constraints.
 * Really just has k (finite stiffness) and lambda (lagrange multipliers)
 */
struct HardConstraintEnergyPool : TombstonePool
{
    Real k_start;           // the initial stiffness
    Real lambda_min;        // lower lambda bound for this constraint (if applicable)
    Real lambda_max;         // upper lambda bound for this constraint (if applicable)
    std::vector<Real> ks;   // the finite stiffnesses of the quadratic energy
    std::vector<Real> lambdas;  // the Lagrange multipliers enforcing the constraints
    std::vector<Real> C_prevs;  // the constraint violation at the end of the previous time step

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
        , ks(capacity, k_start_)
        , lambdas(capacity)
    {}

    /** Add an energy
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy()
    {
        unsigned slot = allocSlot();

        // initialize lambda and k
        lambdas[slot] = 0;
        ks[slot] = k_start;

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