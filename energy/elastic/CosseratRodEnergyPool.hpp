#pragma once

#include "common/TombstonePool.hpp"
#include "energy/elastic/CosseratRodEnergyDefs.hpp"

namespace Energy
{

/** Pool of memory for the Cosserat rod elastic energies. */
struct CosseratRodEnergyPool : TombstonePool
{
    static constexpr int NumParticlesPerEnergy = 2;
    static constexpr EnergyType Type = EnergyType::COSSERAT_ROD;
    static constexpr StaticEnergyType StaticType = StaticEnergyType::COSSERAT_ROD;
    using SolverType = CosseratRodEnergySolver;

    std::vector<CosseratRodEnergyInfo> data;

    explicit CosseratRodEnergyPool(unsigned capacity)
        : TombstonePool(capacity)
        , data(capacity)
    {

    }

    /** Returns a non-owning view into the memory pool. */
    CosseratRodEnergyPoolView view()
    {
        CosseratRodEnergyPoolView v;
        TombstonePool::view(v);
        v.data = data.data();
        return v;
    }

    /** Add an energy
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(
        const Vec2u& indices,
        const Vec6r& stiffness,
        Real rest_length,
        const Vec3r& precurvature
    )
    {
        unsigned slot = allocSlot();

        data[slot].particle_indices = indices;
        data[slot].stiffness = stiffness;
        data[slot].rest_length = rest_length;
        data[slot].precurvature = precurvature;

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