#pragma once

#include "common/TombstonePool.hpp"

namespace Energy
{

struct NeoHookeanEnergyInfo
{
    Vec4u particle_indices;     // indices of particles in the element
    Real lambda;                // first Lame parameter
    Real mu;                    // second Lame parameter
    Real kd;                    // strain-rate damping coefficient
    Mat3r Q;                    // rest state matrices (F = XQ^T)
    Real rest_volume;           // rest-state volume of the elements
    Mat3r E_prev;               // previous Green strain (used for strain-rate damping)
};

/** Pool of memory for the Neo-Hookean energies.
 */
struct NeoHookeanEnergyPool : TombstonePool
{
    static constexpr int NumParticlesPerEnergy = 4; // number of particles per energy
    static constexpr EnergyType Type = EnergyType::NEO_HOOKEAN;   // type of energy in the EnergyType enum
    using Solver = NeoHookeanEnergySolver;  // the solver struct for this pool
     
    std::vector<NeoHookeanEnergyInfo> data;

    explicit NeoHookeanEnergyPool(unsigned capacity)
        : TombstonePool(capacity)
        , data(capacity)
        // , particle_indices(capacity)
        // , lambdas(capacity)
        // , mus(capacity)
        // , kds(capacity)
        // , Qs(capacity)
        // , rest_volumes(capacity)
    {

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
        unsigned slot = allocSlot();

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

    /** Remove an energy
     * @param slot : the index of the energy in the pool to remove
     */
    void removeEnergy(unsigned slot)
    {
        freeSlot(slot);
    }
};

} // namespace Energy