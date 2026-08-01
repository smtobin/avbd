#pragma once

#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

struct GroundCollisionEnergyInfo : HardConstraintEnergyInfo
{
    Vec1u particle_indices;     // particle indices for each constraint
};

struct GroundCollisionEnergyPool : HardConstraintEnergyPool<GroundCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 1; // number of particles per constraint
    static constexpr EnergyType Type = EnergyType::GROUND_COLLISION; // type of energy in the EnergyType enum
    static constexpr StaticEnergyType StaticType = StaticEnergyType::GROUND_COLLISION; // type of energy in the StaticEnergyType enum
    using SolverType = GroundCollisionEnergySolver;     // solver class type

    explicit GroundCollisionEnergyPool(unsigned capacity)
        : HardConstraintEnergyPool(capacity, 1e2, 0)
    {

    }

    /** Add an energy
     * @param particle_index : the index of the particle in the particle pool
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(unsigned particle_index)
    {
        // parent will call allocSlot()
        unsigned slot = HardConstraintEnergyPool::addEnergy();

        // particle_indices[slot][0] = particle_index;
        data[slot].particle_indices[0] = particle_index;

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