#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "energy/CosseratRodEnergyPool.hpp"

namespace Energy
{

/** Implements linear Cosserat rod finite-element energies */
struct CosseratRodEnergySolver
{
    /** Public typedefs */
    using PoolType = CosseratRodEnergyPool;
    static constexpr bool SupportsPositional = false;
    static constexpr bool SupportsOriented = true;

    /** Energy function */
    static Real energy(
        unsigned e_idx,
        const CosseratRodEnergyPool& energies,
        ParticlePool& particles,
        Real /* dt */
    )
    {
        /** TODO: (08/20/26) */
    }

    /** Required - does nothing */
    static void updateAfterIteration(
        unsigned /* e_idx */,
        const CosseratRodEnergyPool& /* energies */,
        ParticlePool& /* particles */
    )
    {

    }

    /** Required - does nothing */
    static void updateAfterTimeStep(
        unsigned e_idx,
        CosseratRodEnergyPool& energies,
        ParticlePool& particles
    )

    /** Computes the Hessian and gradient for a specified particle affected by this energy.
     * Updates the accumulated vertex Hessian and gradients.
     * @param e_idx : the energy index
     * @param energies : the memory pool for the energies
     * @param particles : the memory pool for the particles
     * @param local_idx : the "local" index of the particle in the energy
     * @param particle_H : the current accumulated particle Hessian - this is updated by this function
     * @param particle_G : the current accumulated particle gradient - this is updated by this function
     * @param dt : the time step
    */
    static void accumulate(
        unsigned e_idx,
        const CosseratRodEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat6r& particle_H,
        Mat3r& particle_G,
        Real /* dt */
    )
    {
        /** TODO: (08/20/26) */
    }


};

} // namespace Energy