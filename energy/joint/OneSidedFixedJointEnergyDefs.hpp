#pragma once

#include "energy/HardConstraintEnergyDefs.hpp"

namespace Energy
{

/** Required data for one sided fixed joint energy. */
struct OneSidedFixedJointEnergyInfo : HardConstraintEnergyInfo<6>
{
    Vec1u particle_indices;     // particle indices for each constraint

    Vec3r position_offset;           // positional offset from particle to joint frame
    Quaternion rotation_offset;      // rotational offset from particle to joint frame
    
    Vec3r ref_position;         // reference position
    Quaternion ref_rotation;    // reference orientation
};

///////////////////////////////////////////////////////////////////////////////

/** A non-owning view into the one-sided fixed joint energy storage. */
struct OneSidedFixedJointEnergyPoolView : HardConstraintEnergyPoolView<OneSidedFixedJointEnergyInfo>
{
    using PoolType = OneSidedFixedJointEnergyPool;
};

} // namespace Energy