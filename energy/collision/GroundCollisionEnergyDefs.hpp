#pragma once

#include "energy/collision/CollisionConstraintEnergyDefs.hpp"

namespace Energy
{

/** Required information for ground collision energies. */
struct GroundCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec1u particle_indices;     // particle indices for each constraint
    Real cp_x;
    Real cp_z;
};

//////////////////////////////////////////////////////////////////////////////////

/** A non-owning view into ground collision energy storage */
struct GroundCollisionEnergyPoolView : CollisionConstraintEnergyPoolView<GroundCollisionEnergyInfo>
{
};

} // namespace Energy