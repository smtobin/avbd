#pragma once

#include "CollisionConstraintEnergyPool.hpp"

namespace Energy
{

template <typename CollisionEnergyPool, typename CollisionConstraintSolver>
struct CollisionConstraintEnergySolver
{
    /** Updates the stiffness and Lagrange multiplier after the iteration
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     * @param particles : the simulation particle memory pool
     */
    static void updateAfterIteration(
        unsigned c_idx,
        CollisionEnergyPool& energies,
        ParticlePool& particles
    )
    {
        // evaluate the constraint
        Vec3r C = CollisionConstraintSolver::evaluateConstraint(c_idx, energies, particles);

        // subtract off previous constraint violation
        Vec3r C_corr = C - CONSTRAINT_ALPHA * energies.data[c_idx].C_prev;
        
        // extract the current stiffness and Lagrange multiplier
        Vec3r& k = energies.data[c_idx].k;
        Vec3r& lambda = energies.data[c_idx].lambda;

        Vec3r lambda_p = k.cwiseProduct(C_corr) + lambda;

        // update normal stiffness - equation (12)
        // for collision constraints, lambda min = 0
        if (lambda_p[2] > 0)
        {
            // additive stiffness update
            k[2] += STIFFNESS_BETA * C_corr[2];

            // store the current constraint violation
            energies.data[c_idx].C_prev[2] = C[2];
        }

        // update normal lambda - equation (11)
        lambda[2] = std::max(Real(0), lambda_p[2]);

        // update tangent and binromal stiffness - equation (12)
        Vec2r lambda_tb_p(lambda_p[0], lambda_p[1]);
        Real lambda_tb_p_mag = std::sqrt(lambda_p[0]*lambda_p[0] + lambda_p[1]*lambda_p[1]);
        Real mu = 0.2;
        /** TODO: (08/04/26) make coeff of friction a part of data.
         * TODO: (08/04/26) static + dynamic friction?
         */
        Real lambda_tb_max = mu*lambda_p[2];
        if (lambda_tb_p_mag < lambda_tb_max)
        {
            k[0] += STIFFNESS_BETA * C_corr[0];
            k[1] += STIFFNESS_BETA * C_corr[1];

            energies.data[c_idx].C_prev[0] = C[0];
            energies.data[c_idx].C_prev[1] = C[1];
        }
        else
        {
            Vec2r lambda_tb_clamped = lambda_tb_p / (lambda_tb_p_mag + 1e-8) * lambda_tb_max;
            lambda[0] = lambda_tb_clamped[0];
            lambda[1] = lambda_tb_clamped[1];
        }

        
    }

    /** Updates the stiffness and Lagrange multiplier after a full time step.
     * Implements equation (19) from the AVBD paper.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     */
    static void updateAfterTimeStep(
        unsigned c_idx,
        CollisionEnergyPool& energies,
        ParticlePool& /* particles */
    )
    {
        energies.data[c_idx].lambda = CONSTRAINT_ALPHA * STIFFNESS_GAMMA * energies.data[c_idx].lambda;
        energies.data[c_idx].k = (STIFFNESS_GAMMA * energies.data[c_idx].k).cwiseMax(energies.k_start);
    }

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
        const CollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat3r& particle_H,
        Vec3r& particle_G,
        Real dt
    )
    {
        // evaluate the constraint, gradients, and Hessians for the constraint for each particle involved
        Vec3r C_grad_n, C_grad_t, C_grad_b;
        Mat3r C_hess_n, C_hess_t, C_hess_b;
        Vec3r C_raw;
        CollisionConstraintSolver::constraintGradientHessian(
            e_idx, energies, particles, local_idx, // inputs
            C_raw, 
            C_grad_n, C_grad_t, C_grad_b,
            C_hess_n, C_hess_t, C_hess_b  // outputs
        );
        // alpha correction on the constraint violation
        Vec3r C = C_raw - CONSTRAINT_ALPHA * energies.data[e_idx].C_prev;

        const Vec3r& k = energies.data[e_idx].k;
        const Vec3r& lambda = energies.data[e_idx].lambda;

        // Lagrange multiplier for normal
        Real lambda_p_n = std::max(Real(0), k[2]*C[2] + lambda[2]);

        // Lagrange multipliers for tangent and binormal
        Vec2r lambda_p_tb(k[0]*C[0] + lambda[0], k[1]*C[1] + lambda[1]);
        // clamp magnitude
        Real lambda_p_tb_mag = lambda_p_tb.norm();
        Real mu = 0.2;
        Real lambda_p_tb_max = mu*lambda_p_n;
        if (lambda_p_tb_mag > lambda_p_tb_max)
        {
            lambda_p_tb = lambda_p_tb / (lambda_p_tb.norm() + 1e-8) * lambda_p_tb_max;
        }
        
        
        /** TODO: (08/04/26) After we've clamped the lambdas, these will always be false...? */
        // stiffness rescaling for normal - equation (14)
        Real k_scaled_n = k[2];
        if (lambda_p_n < 0 && std::abs(C[2]) > 1e-12)
            k_scaled_n =  -lambda[2] / C[2];
        
        // stiffness rescling for tangent and binormal - equation (14)
        Vec2r k_scaled_tb(k[0], k[1]);
        Vec2r C_tb(C[0], C[1]);
        if (lambda_p_tb.norm() > lambda_p_tb_max && C_tb.norm() > 1e-12)
        {
            k_scaled_tb = k_scaled_tb / k_scaled_tb.norm() * (lambda_p_tb_max - lambda_p_tb.norm()) / C_tb.norm();
        }
        // gradient
        Vec3r grad_n = lambda_p_n * C_grad_n;
        Vec3r grad_t = lambda_p_tb[0] * C_grad_t;
        Vec3r grad_b = lambda_p_tb[1] * C_grad_b;

            // Hessian
            /** TODO: do diagonalization of Hessian component
             * 
             * 
             * 
             * 
             */
            

        Mat3r hess_n = (k_scaled_n * C[2] + lambda[2]) * C_hess_n +
            k_scaled_n * C_grad_n* C_grad_n.transpose();
        Mat3r hess_t = (k_scaled_tb[0] * C[0] + lambda[0]) * C_hess_t +
            k_scaled_tb[0] * C_grad_t * C_grad_t.transpose();
        Mat3r hess_b = (k_scaled_tb[1] * C[1] + lambda[1]) * C_hess_b +
            k_scaled_tb[1] * C_grad_b * C_grad_b.transpose();

        particle_H += hess_n + hess_t + hess_b;
        particle_G += grad_n + grad_t + grad_b;
    }
};

} // namespace Energy