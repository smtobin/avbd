#pragma once

#include "common/TombstonePool.hpp"

namespace Energy
{

/** Pool of memory for the Neo-Hookean energies.
 */
struct NeoHookeanEnergyPool : TombstonePool
{
    static constexpr int NumParticlesPerEnergy = 4; // number of particles per energy
    static constexpr EnergyType Type = EnergyType::NEO_HOOKEAN;   // type of energy in the EnergyType enum

    std::vector<Vec4u> particle_indices;    // indices of particles in the element
    std::vector<Real> lambdas;              // first Lame parameter
    std::vector<Real> mus;                  // second Lame parameter
    std::vector<Real> kds;                  // strain-rate damping coefficient
    std::vector<Mat3r> Qs;                  // rest state matrices (F = XQ^T)
    std::vector<Real> rest_volumes;         // rest-state volume of the elements

    explicit NeoHookeanEnergyPool(unsigned capacity)
        : TombstonePool(capacity)
        , particle_indices(capacity)
        , lambdas(capacity)
        , mus(capacity)
        , kds(capacity)
        , Qs(capacity)
        , rest_volumes(capacity)
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

        // copy over the input information
        particle_indices[slot] = indices;
        lambdas[slot] = lambda;
        mus[slot] = mu;
        kds[slot] = kd;
        Qs[slot] = Q;
        rest_volumes[slot] = rest_volume;

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