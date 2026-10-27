#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"

namespace Energy
{

/** Required data for all hard constraints (other than collisions). */
template <int N>
struct HardConstraintEnergyInfo
{
    VecNr<N> k;         // the finite stiffnesses of the quadratic energy
    VecNr<N> lambda;    // the Lagrange multipliers enforcing the constraints
    VecNr<N> C_prev;    // the constraint violation at the end of the previous time step
};

///////////////////////////////////////////////////////////////////////////////////////////////////

/** A non-owning view into the hard constraint energy storage. */
template <typename EnergyInfo>
struct HardConstraintEnergyPoolView : TombstonePoolView
{
    EnergyInfo* data;
};

} // namespace Energy