#pragma once

#include "common/TombstonePool.hpp"

namespace Energy
{

template <int N>
struct HardConstraintEnergyInfo
{
    VecNr<N> k;         // the finite stiffnesses of the quadratic energy
    VecNr<N> lambda;    // the Lagrange multipliers enforcing the constraints
    VecNr<N> C_prev;    // the constraint violation at the end of the previous time step
};

/** Pool of memory for HardConstraintEnergies.
 * Base class for all energies associated with hard constraints.
 * Really just has k (finite stiffness) and lambda (lagrange multipliers)
 */
template <typename EnergyInfo, int N>
struct HardConstraintEnergyPool : TombstonePool
{
    static constexpr int Dim = N;
    
    VecNr<N> k_start;           // the initial stiffness
    VecNr<N> lambda_min;        // lower lambda bound for this constraint (if applicable)
    VecNr<N> lambda_max;         // upper lambda bound for this constraint (if applicable)
    std::vector<EnergyInfo> data;

    explicit HardConstraintEnergyPool(
        unsigned capacity,
        VecNr<N> k_start_,
        VecNr<N> lambda_min_= VecNr<N>::Constant(std::numeric_limits<Real>::lowest()),
        VecNr<N> lambda_max_= VecNr<N>::Constant(std::numeric_limits<Real>::max())
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
        data[slot].lambda = VecNr<N>::Zero();
        data[slot].k = k_start;
        data[slot].C_prev = VecNr<N>::Zero();

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