#pragma once

#include "common/TombstonePool.hpp"

namespace Energy
{

struct CosseratRodEnergyInfo
{
    Vec2u particle_indices;     // the indices of the oriented particles in the rod element
    Vec6r stiffness;            // the diagonal of the body-frame stiffness matrix (GA, GA, EA, EI, EI, GJ)
    Real rest_length;           // initial length of the element
    Vec3r precurvature;         // initial curvature of the element
};

/** Pool of memory for the Cosserat rod elastic energies. */
struct CosseratRodEnergyPool : TombstonePool
{
    static constexpr int NumParticlePerEnergy = 2;
    static constepxr EnergyType = EnergyType::COSSERAT_ROD;
    static constexpr StaticEnergyType StaticType = StaticEnergyType::COSSERAT_ROD;
    using SolverType = CosseratRodEnergySolver;

    std::vector<CosseratRodEnergyInfo> data;

    explicit CosseratRodEnergyPool(unsigned capacity)
        : TombstonePool(capacity)
        , data(capacity)
    {

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