#pragma once

#include "common/common.hpp"
#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

template <typename EnergyPool, typename ConstraintSolver>
struct HardConstraintEnergySolver
{
    /** Updates the stiffness and Lagrange multiplier after the iteration
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     * @param particles : the simulation particle memory pool
     */
    static void updateAfterIteration(
        unsigned c_idx,
        const EnergyPool& energies,
        ParticlePool& particles
    )
    {
        // evaluate the constraint
        Real C = ConstraintSolver::evaluateConstraint(c_idx, energies, particles);

        // subtract off previous constraint violation
        Real C_corr = C - CONSTRAINT_ALPHA * energies.C_prevs[c_idx];
        
        // extract the current stiffness and Lagrange multiplier
        Real k = energies.ks[c_idx];
        Real lambda = energies.lambdas[c_idx];

        Real lambda_p = k * C_corr + lambda;

        // update stiffness - equation (12)
        if (lambda_p > energies.lambda_min && lambda_p < energies.lambda_max)
        {
            // additive stiffness update
            energies.ks[c_idx] += STIFFNESS_BETA * C_corr;

            // store the current constraint violation
            energies.C_prevs[c_idx] = C;
        }

        // update lambda - equation (11)
        energies.lambdas[c_idx] = std::max(lambda_min, std::min(lambda_max, lambda_p));
    }

    /** Updates the stiffness and Lagrange multiplier after a full time step.
     * Implements equation (19) from the AVBD paper.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     */
    static void updateAfterTimeStep(
        unsigned c_idx,
        const EnergyPool& energies
    )
    {
        energies.lambdas[c_idx] = CONSTRAINT_ALPHA * STIFFNESS_GAMMA * energies.lambdas[c_idx];
        energies.ks[c_idx] = std::max(STIFFNESS_GAMMA * energies.ks[c_idx], energies.k_start);
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
        const EnergyPool& energies,
        ParticlePool& particles,
        Mat3r* vertex_Hs,
        Vec3r* vertex_Gs,
        Real dt
    )
    {
        // evaluate the constraint, gradients, and Hessians for the constraint for each particle involved
        Vec3r grads[EnergyPool::NumParticlesPerConstraint];
        Mat3r hessians[EnergyPool::NumParticlesPerConstraint];
        Real C_raw;
        ConstraintSolver::constraintGradientHessian(
            c_idx, energies, particles, // inputs
            C_raw, grads, hessians  // outputs
        );
        // alpha correction on the constraint violation
        Real C = C_raw - CONSTRAINT_ALPHA * energies.C_prevs[c_idx];

        Real k = energies.ks[c_idx];
        Real lambda = energies.lambdas[c_idx];

        Real lambda_p = std::max(energies.lambda_min, std::min(energies.lambda_max, k*C + lambda));
        
        // stiffness rescaling - equation (14)
        Real k_scaled = k;
        if (lambda_p < energies.lambda_min && std::abs(C) > 1e-12)
            k_scaled = (energies.lambda_min - lambda) / C;
        else if (lambda_p > energies.lambda_max && std::abs(C) > 1e-12)
            k_scaled = (energies.lambda_max - lambda) / C;
        
        const auto& indices = energies.particles;
        for (int i = 0; i < EnergyPool::NumParticlesPerConstraint; i++)
        {
            // gradient
            Vec3r grad_i = lambda_p * grads[i];

            // Hessian
            /** TODO: do diagonalization of Hessian component
             * 
             * 
             * 
             * 
             */
            

            Mat3r hess_i = (k_scaled * C + lambda) * hessians[i] +
                k_scaled * grads[i] * grads[i].transpose();
            
            vertex_Hs[indices[i]] += hess_i;
            vertex_Gs[indices[i]] += grad_i;
        }
        
    }
};

} // namespace Energy