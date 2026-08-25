#pragma once

#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

struct OneSidedFixedJointEnergyInfo : HardConstraintEnergyInfo<6>
{
    Vec1u particle_indices;     // particle indices for each constraint

    Vec3r position_offset;           // positional offset from particle to joint frame
    Quaternion rotation_offset;      // rotational offset from particle to joint frame
    
    Vec3r ref_position;         // reference position
    Quaternion ref_rotation;    // reference orientation
};

struct OneSidedFixedJointEnergyPool : HardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>
{
    static constexpr int NumParticlesPerEnergy = 1; // number of particles per constraint
    static constexpr EnergyType Type = EnergyType::ONE_SIDED_FIXED_JOINT; // type of energy in the EnergyType enum
    static constexpr StaticEnergyType StaticType = StaticEnergyType::ONE_SIDED_FIXED_JOINT; // type of energy in the StaticEnergyType enum
    using SolverType = OneSidedFixedJointEnergySolver;     // solver class type

    explicit OneSidedFixedJointEnergyPool(unsigned capacity)
        : HardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>(capacity, Vec6r::Constant(1))
    {

    }

    /** Add an energy
     * @param particle_index : the index of the particle in the particle pool
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(
        unsigned particle_index, 
        const Vec3r& pos_offset,
        const Quaternion& rot_offset,
        const Vec3r& ref_position,
        const Quaternion& ref_rotation
    )
    {
        // parent will call allocSlot()
        unsigned slot = HardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>::addEnergy();

        // particle_indices[slot][0] = particle_index;
        data[slot].particle_indices[0] = particle_index;
        data[slot].position_offset = pos_offset;
        data[slot].rotation_offset = rot_offset;
        data[slot].ref_position = ref_position;
        data[slot].ref_rotation = ref_rotation;

        return slot;
    }
};

} // namespace Energy