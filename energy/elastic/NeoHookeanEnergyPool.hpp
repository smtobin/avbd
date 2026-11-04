#pragma once

#include "common/PoolBase.hpp"
#include "energy/elastic/NeoHookeanEnergyDefs.hpp"

namespace Energy
{

/** Pool of memory for the Neo-Hookean energies.
 */
struct NeoHookeanEnergyPool : PoolBase
{
    static constexpr int NumParticlesPerEnergy = 4; // number of particles per energy
    static constexpr EnergyType Type = EnergyType::NEO_HOOKEAN;   // type of energy in the EnergyType enum
    static constexpr StaticEnergyType StaticType = StaticEnergyType::NEO_HOOKEAN; // type of energy in the StaticEnergyType enum
    using SolverType = NeoHookeanEnergySolver;  // the solver struct for this pool
    using ViewType = NeoHookeanEnergyPoolView;  // the view struct for this pool
     
    std::vector<NeoHookeanEnergyInfo> data;

    explicit NeoHookeanEnergyPool(unsigned capacity)
        : PoolBase(capacity)
        , data(capacity)
    {

    }

    /** Returns a non-owning view into the memory pool. */
    NeoHookeanEnergyPoolView view()
    {
        NeoHookeanEnergyPoolView v;
        PoolBase::view(v);
        v.data = data.data();
        return v;
    }

    /** Add an energy
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(
        const Vec4u& indices,
        Real lambda,
        Real mu,
        Real kd,
        const Mat3r& Q,
        Real rest_volume
    )
    {
        unsigned slot = allocate();

        data[slot].particle_indices = indices;
        data[slot].lambda = lambda;
        data[slot].mu = mu;
        data[slot].kd = kd;
        data[slot].Q = Q;
        data[slot].rest_volume = rest_volume;

        // copy over the input information
        // particle_indices[slot] = indices;
        // lambdas[slot] = lambda;
        // mus[slot] = mu;
        // kds[slot] = kd;
        // Qs[slot] = Q;
        // rest_volumes[slot] = rest_volume;

        return slot;
    }
};

} // namespace Energy