#pragma once

#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

struct GroundCollisionEnergyPool : HardConstraintEnergyPool
{
    static int NumParticlesPerConstraint = 1; // number of particles per constraint
    static EnergyType Type = EnergyType::GROUND_COLLISION; // type of energy in the EnergyType enum

    std::vector<Eigen::Vector<unsigned, NumParticlesPerConstraint>> particles;     // particle indices for each constraint

    explicit GroundCollisionEnergyPool(unsigned capacity)
        : HardConstraintEnergyPool(capacity)
        , particles(capacity)
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

        particles[slot] = particle_index;

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