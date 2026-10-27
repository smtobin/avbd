#pragma once

#include "common/TombstonePool.hpp"
#include "energy/HardConstraintEnergyDefs.hpp"

namespace Energy
{

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
        VecNr<N> lambda_min_= VecNr<N>(std::numeric_limits<Real>::lowest()),
        VecNr<N> lambda_max_= VecNr<N>(std::numeric_limits<Real>::max())
    )
        : TombstonePool(capacity)
        , k_start(k_start_)
        , lambda_min(lambda_min_)
        , lambda_max(lambda_max_)
        , data(capacity)
    {}

    /** Returns a non-owning view into the memory pool. */
    void view(HardConstraintEnergyPoolView<EnergyInfo, N>& v)
    {
        TombstonePool::view(v);
        v.data = data.data();
        v.k_start = k_start;
        v.lambda_min = lambda_min;
        v.lambda_max = lambda_max;
    }

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