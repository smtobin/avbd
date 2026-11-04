#pragma once

#include "common/PoolBase.hpp"
#include "energy/elastic/CosseratRodEnergyDefs.hpp"

namespace Energy
{

/** Pool of memory for the Cosserat rod elastic energies. */
struct CosseratRodEnergyPool : PoolBase
{
    static constexpr int NumParticlesPerEnergy = 2;
    static constexpr EnergyType Type = EnergyType::COSSERAT_ROD;
    static constexpr StaticEnergyType StaticType = StaticEnergyType::COSSERAT_ROD;
    using SolverType = CosseratRodEnergySolver;
    using ViewType = CosseratRodEnergyPoolView;

    std::vector<CosseratRodEnergyInfo> data;

    explicit CosseratRodEnergyPool(unsigned capacity)
        : PoolBase(capacity)
        , data(capacity)
    {

    }

    /** Returns a non-owning view into the memory pool. */
    CosseratRodEnergyPoolView view()
    {
        CosseratRodEnergyPoolView v;
        PoolBase::view(v);
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
        unsigned slot = allocate();

        data[slot].particle_indices = indices;
        data[slot].stiffness = stiffness;
        data[slot].rest_length = rest_length;
        data[slot].precurvature = precurvature;

        return slot;
    }
};



} // namespace Energy