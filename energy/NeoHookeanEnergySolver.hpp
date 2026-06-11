#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "energy/NeoHookeanEnergyPool.hpp"

namespace Energy
{

/** Implements the stable Neo-Hookean per-element energies, seen in Macklin et al. 2021 */
struct NeoHookeanEnergySolver
{
    static void computeF(
        const Vec3r& v1,
        const Vec3r& v2,
        const Vec3r& v3, 
        const Vec3r& v4,
        const Mat3r& Q)
    {
        Mat3r X;
        X.col(0) = v1 - v4;
        X.col(1) = v2 - v4;
        X.col(2) = v3 - v4;

        return X*Q;
    }

    /** Computes the Hessian and gradient for all of the particles affected by this constraint.
     * Updates the accumulated vertex Hessian and gradients.
     * @param c_idx : the energy index
     * @param energies : the memory pool for the energies
     * @param particles : the memory pool for the particles
     * @param vertex_Hs : a memory pool storing the per-vertex Hessians - this is updated by this function
     * @param vertex_Gs : a memory pool storing the per-vertex gradients - this is updated by this function
     * @param dt : the time step
     */
    static void accumulateHessianGradient(
        unsigned c_idx,
        const NeoHookeanEnergyPool& energies,
        ParticlePool& particles,
        Mat3r* vertex_Hs,
        Vec3r* vertex_Gs,
        Real dt
    )
    {
        // extract data from energy pool
        const Vec4i& indices = energies.particle_indices[c_idx];
        Real V = energies.rest_volumes[c_idx];
        const Mat3r& Q = energies.Qs[c_idx];
        Real lambda = energies.lambdas[c_idx];
        Real mu = energies.mus[c_idx];
        Real kd = energies.kds[c_idx];

        // compute F for this timestep and the previous timestep
        Mat3r F = computeF(
            particles.positions[indices[0]],
            particles.positions[indices[1]],
            particles.positions[indices[2]],
            particles.positions[indices[3]],
            Q);
        Mat3r F_prev = computeF(
            particles.previous_positions[indices[0]],
            particles.previous_positions[indices[1]],
            particles.previous_positions[indices[2]],
            particles.previous_positions[indices[3]],
            Q);
        
        // compute rate of change of Green strain
        Mat3r E = F.transpose() * F - Mat3r::Identity();
        Mat3r E_prev = F_prev.transpose() * F_prev - Mat3r::Identity();
        Mat3r E_dot = 1/dt * (E - E_prev);

        // hydrostatic gradient
        Mat3r F_cross;
        F_cross.col(0) = F.col(1).cross(F.col(2));
        F_cross.col(1) = F.col(2).cross(F.col(0));
        F_cross.col(2) = F.col(0).cross(F.col(1));

        Mat3r hyd_grad_full = V*lambda * (F.determinant() - mu/lambda - 1) * F_cross * Q.transpose();

        // deviatoric gradient
        Mat3r dev_grad_full = V*mu * F * Q.transpose();

        // compute the Hessian and gradient for the first 3 particles
        for (int i = 0; i < 3; i++)
        {
            Vec3r qi = Q.row(i);
            Vec3r Fqi = F*qi;

            // gradient
            Vec3r hyd_grad_i = hyd_grad_full.col(i);
            Vec3r dev_grad_i = dev_grad_full.col(i);
            Vec3r damp_grad_i = 2/dt * V * kd * F * E_dot * qi;
            vertex_Gs[indices[i]] += hyd_grad_i + dev_grad_i + damp_grad_i;

            // Hessian
            Mat3r hyd_hess_i = V * lambda * detF_grad_full.col(i) * detF_grad_full.col(i).transpose();
            Mat3r dev_hess_i = V * mu * qi.squaredNorm() * Mat3r::Identity();
            Mat3r damp_hess_i = 2/dt * V * kd * ( (qi.transpose() * E_dot * qi) * Mat3r::Identity() + 1/dt * (Fqi * Fqi.transpose() + F*F.transpose() * qi.squaredNorm()) );
            vertex_Hs[indices[i]] += hyd_hess_i + dev_hess_i + damp_hess_i;
        }

        // compute Hessian and gradient for the 4th particle
        // gradient
        Vec3r q4 = Q.row(0) + Q.row(1) + Q.row(2);
        Vec3r hyd_grad_4 = -hyd_grad_full.col(0) - hyd_grad_full.col(1) - hyd_grad_full.col(2);
        Vec3r dev_grad_4 = -dev_grad_full.col(0) - dev_grad_full.col(1) - dev_grad_full.col(2);
        Vec3r damp_grad_4 = -2/dt * V * kd * F * E_dot * q4;
        vertex_Gs[indices[i]] += hyd_grad_4 + dev_grad_4 + damp_grad_4;

        // Hessian
        Vec3r a3 =
            -detF_grad_full.col(0)
            -detF_grad_full.col(1)
            -detF_grad_full.col(2);

        Mat3r hyd_hess =
            V*lambda * a3 * a3.transpose();

        Mat3r dev_hess =
            V*mu * q4.squaredNorm() * Mat3r::Identity();

        Vec3r Fq4 = F * q4;
        Mat3r damp_hess = 2/dt * V * kd * ( (q4.transpose() * E_dot * q4) * Mat3r::Identity() + 1/dt * (Fq4 * Fq4.transpose() + F*F.transpose() * q4.squaredNorm()) );
        vertex_Hs[indices[i]] += hyd_hess + dev_hess + damp_hess;
    }
};

} // namespace Energy