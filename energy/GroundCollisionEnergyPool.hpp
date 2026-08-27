#pragma once

#include "energy/CollisionConstraintEnergyPool.hpp"

namespace Energy
{

struct GroundCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec1u particle_indices;     // particle indices for each constraint
    Real cp_x;
    Real cp_z;
};

struct GroundCollisionEnergyPool : CollisionConstraintEnergyPool<GroundCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 1; // number of particles per constraint
    static constexpr EnergyType Type = EnergyType::GROUND_COLLISION; // type of energy in the EnergyType enum
    static constexpr StaticEnergyType StaticType = StaticEnergyType::GROUND_COLLISION; // type of energy in the StaticEnergyType enum
    using SolverType = GroundCollisionEnergySolver;     // solver class type

    explicit GroundCollisionEnergyPool(unsigned capacity)
        : CollisionConstraintEnergyPool(capacity)
    {

    }

    /** Add an energy
     * @param particle_index : the index of the particle in the particle pool
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(unsigned particle_index, Real k_start, Real mu_s, Real mu_k)
    {
        // parent will call allocSlot()
        unsigned slot = CollisionConstraintEnergyPool::addEnergy(k_start, mu_s, mu_k);

        // particle_indices[slot][0] = particle_index;
        data[slot].particle_indices[0] = particle_index;
        data[slot].cp_x = 0;
        data[slot].cp_z = 0;

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