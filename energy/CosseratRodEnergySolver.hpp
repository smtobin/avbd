#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "common/Math.hpp"
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
        strain.block<3,1>(0,0) = 1/rest_length * (q_mid.conjugate() * (p2 - p1)) - Vec3r(0,0,1);

        // bending strain
        strain.block<3,1>(3,0) = 1.0/rest_length * q_diff - precurvature;

        return 0.5 * strain.transpose() * stiffness.asDiagonal() * strain;
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
    {}

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
        Vec6r& particle_G,
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
        strain.block<3,1>(0,0) = 1/rest_length * (q_mid.conjugate() * (p2 - p1)) - Vec3r(0,0,1);

        // bending strain
        strain.block<3,1>(3,0) = 1.0/rest_length * q_diff - precurvature;

        Real inv_length = 1.0/rest_length;

        Mat3r gam_inv = Math::ExpMap_InvRightJacobian(q_diff);

        // theta(s)
        Vec3r theta = s_hat * q_diff;

        // precompute useful quantities for gradients of shear strain
        Mat3r exp_theta = Math::Exp_so3(theta);
        Mat3r gam_theta = Math::ExpMap_RightJacobian(theta);
        Mat3r R = q_mid.toRotationMatrix();


        Vec3r dp_ds = inv_length * (p2 - p1);
        Mat3r RT_dp_ds_R = Math::Skew3(R.transpose() * dp_ds);

        Mat3r dtheta_dRi, dtheta_ds_dRi;
        Mat6r strain_grad;
        strain_grad.block<3,3>(3,0) = Mat3r::Zero();
        if (local_idx == 0)
        {
            dtheta_dRi = -s_hat * gam_inv.transpose();
            dtheta_ds_dRi = -inv_length * gam_inv.transpose();

            strain_grad.block<3,3>(0,0) = -inv_length * R.transpose();
            strain_grad.block<3,3>(0,3) =  RT_dp_ds_R * (exp_theta.transpose() + gam_theta * dtheta_dRi);
            strain_grad.block<3,3>(3,3) = dtheta_ds_dRi;
        }
        else
        {
            dtheta_dRi = s_hat * gam_inv;
            dtheta_ds_dRi = inv_length * gam_inv;

            strain_grad.block<3,3>(0,0) = inv_length * R.transpose();
            strain_grad.block<3,3>(0,3) = RT_dp_ds_R * gam_theta * dtheta_dRi;
            strain_grad.block<3,3>(3,3) = dtheta_ds_dRi;
        }

        particle_G += strain.transpose() * (stiffness.asDiagonal() * strain_grad);
        /** TODO: (08/23/26) Is Gauss-Newton approximation good enough? */
        particle_H += strain_grad.transpose() * stiffness.asDiagonal() * strain_grad; // Gauss-Newton approximation - ignoring constraint Hessian term
    }


};

} // namespace Energy