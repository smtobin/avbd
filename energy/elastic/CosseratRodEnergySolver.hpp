#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/Math.hpp"
#include "energy/elastic/CosseratRodEnergyPool.hpp"
#include "energy/elastic/CosseratRodEnergyHelpers.hpp"

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
        const CosseratRodEnergyPoolView& energies,
        ParticlePoolView& particles,
        Real /* dt */
    )
    {
        const CosseratRodEnergyInfo& info = energies.data[e_idx];
        const Vec2u& indices = info.particle_indices;
        const Quaternion& q1 = particles.rotation(indices[0]);
        const Quaternion& q2 = particles.rotation(indices[1]);
        const Vec3r& p1 = particles.positions[indices[0]];
        const Vec3r& p2 = particles.positions[indices[1]];
        Real rest_length = info.rest_length;
        const Vec3r& precurvature = info.precurvature;
        const Vec6r& stiffness = info.stiffness;
        Real s_hat = 0.5;

        Vec3r q_diff = Math::Minus_S3(q2, q1);
        Vec6r strain;

        // shear strain
        Quaternion q_mid = Math::Plus_S3(q1, s_hat * q_diff);
        strain.head<3>() = 1/rest_length * (q_mid.conjugate() * (p2 - p1)) - Vec3r(0,0,1);

        // bending strain
        strain.tail<3>() = 1.0/rest_length * q_diff - precurvature;

        return 0.5 * rest_length * strain.dot(stiffness.cwiseProduct(strain));
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
        unsigned /* e_idx */,
        CosseratRodEnergyPoolView& /* energies */,
        ParticlePoolView& /* particles */
    )
    {}

    /** Computes the Hessian and gradient for a specified particle affected by this energy.
     * Updates the accumulated vertex Hessian and gradients.
     * @param e_idx : the energy index
     * @param energies : the memory pool for the energies
     * @param particles : the memory pool for the particles
     * @param local_idx : the "local" index of the particle in the energy
     * @param H : (output) the Hessian for this energy
     * @param G : (output) the gradient for this energy
     * @param dt : the time step
    */
    static void accumulate(
        unsigned e_idx,
        const CosseratRodEnergyPoolView& energies,
        ParticlePoolView& particles,
        unsigned local_idx,
        Mat6r& H,
        Vec6r& G,
        Real /* dt */
    )
    {
        const CosseratRodEnergyInfo& info = energies.data[e_idx];
        const Vec2u& indices = info.particle_indices;
        const Quaternion& q1 = particles.rotation(indices[0]);
        const Quaternion& q2 = particles.rotation(indices[1]);
        const Vec3r& p1 = particles.positions[indices[0]];
        const Vec3r& p2 = particles.positions[indices[1]];
        Real rest_length = info.rest_length;
        const Vec3r& precurvature = info.precurvature;
        const Vec6r& stiffness = info.stiffness;

        cosseratRodEnergyAccumulate(
            p1,
            p2,
            q1,
            q2,
            rest_length,
            precurvature,
            stiffness,
            local_idx,
            H,
            G
        );
    }


};

} // namespace Energy