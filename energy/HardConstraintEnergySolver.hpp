#pragma once

#include "common/common.hpp"
#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

template <typename EnergyPool, typename ConstraintSolver>
struct HardConstraintEnergySolver
{
    using Vec = VecNr<EnergyPool::Dim>;

    /** Updates the stiffness and Lagrange multiplier after the iteration
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     * @param particles : the simulation particle memory pool
     */
    static void updateAfterIteration(
        unsigned c_idx,
        EnergyPool& energies,
        ParticlePool& particles
    )
    {
        // evaluate the constraint
        Vec C = ConstraintSolver::evaluateConstraint(c_idx, energies, particles);

        // subtract off previous constraint violation
        Vec C_corr = C - CONSTRAINT_ALPHA * energies.data[c_idx].C_prev;
        
        // extract the current stiffness and Lagrange multiplier
        Vec k = energies.data[c_idx].k;
        Vec lambda = energies.data[c_idx].lambda;

        Vec lambda_p = k.cwiseProduct(C_corr) + lambda;

        // update stiffness - equation (12)
        for (unsigned ci = 0; ci < EnergyPool::Dim; ci++)
        {
            if (lambda_p[ci] > energies.lambda_min[ci] && lambda_p[ci] < energies.lambda_max[ci])
            {
                // additive stiffness update
                energies.data[c_idx].k[ci] += STIFFNESS_BETA * C_corr[ci];
            }

            // update lambda - equation (11)
            energies.data[c_idx].lambda[ci] = std::max(energies.lambda_min[ci], std::min(energies.lambda_max[ci], lambda_p[ci]));
        }
    }

    /** Updates the stiffness and Lagrange multiplier after a full time step.
     * Implements equation (19) from the AVBD paper.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     */
    static void updateAfterTimeStep(
        unsigned c_idx,
        EnergyPool& energies,
        ParticlePool& particles
    )
    {
        energies.data[c_idx].lambda = CONSTRAINT_ALPHA * STIFFNESS_GAMMA * energies.data[c_idx].lambda;
        energies.data[c_idx].k = energies.k_start.cwiseMax(STIFFNESS_GAMMA * energies.data[c_idx].k);
        energies.data[c_idx].C_prev = ConstraintSolver::evaluateConstraint(c_idx, energies, particles);
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
    template <int DOF>
    static void accumulate(
        unsigned e_idx,
        const EnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat3r_or_Mat6r<DOF>& particle_H,
        Vec3r_or_Vec6r<DOF>& particle_G,
        Real dt
    )
    {
        std::cout << "HardconstraintEnergy accumulate! DOF=" << DOF << std::endl;
        // evaluate the constraint, gradients, and Hessians for the constraint for each particle involved
        std::array<Vec3r_or_Vec6r<DOF>, EnergyPool::Dim> C_grad;
        std::array<Mat3r_or_Mat6r<DOF>, EnergyPool::Dim> C_hess;       // constraint Hessian is 3rd order tensor
        Vec C_raw;
        ConstraintSolver::template constraintGradientHessian<DOF>(
            e_idx, energies, particles, local_idx, // inputs
            C_raw, C_grad, C_hess  // outputs
        );
        // alpha correction on the constraint violation
        const Vec C = C_raw - CONSTRAINT_ALPHA * energies.data[e_idx].C_prev;

        const Vec& k = energies.data[e_idx].k;
        const Vec& lambda = energies.data[e_idx].lambda;

        Vec lambda_p_unclamped = k.cwiseProduct(C) + lambda;
        Vec lambda_p = energies.lambda_min.cwiseMax( energies.lambda_max.cwiseMin( lambda_p_unclamped ));
        
        // stiffness rescaling - equation (14)
        Vec k_scaled = k;
        for (unsigned ci = 0; ci < EnergyPool::Dim; ci++)
        {
            if (lambda_p_unclamped[ci] < energies.lambda_min[ci] && std::abs(C[ci]) > 1e-12)
                k_scaled[ci] = std::abs( (energies.lambda_min[ci] - lambda[ci]) / C[ci] );
            else if (lambda_p_unclamped[ci] > energies.lambda_max[ci] && std::abs(C[ci]) > 1e-12)
                k_scaled[ci] = std::abs( (energies.lambda_max[ci] - lambda[ci]) / C[ci] );

            // gradient
            particle_G += C_grad[ci] * lambda_p[ci];

                // Hessian
                /** TODO: do diagonalization of Hessian component
                 * 
                 * 
                 * 
                 * 
                 */
                

            Mat3r_or_Mat6r<DOF> hess = (k_scaled[ci] * C[ci] + lambda[ci]) * C_hess[ci] +
                k_scaled[ci] * C_grad[ci] * C_grad[ci].transpose();

            particle_H += hess;
        }
    }

    static void accumulate(
        unsigned e_idx,
        const EnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat3r& particle_H,
        Vec3r& particle_G,
        Real dt
    )
    {
        accumulate<3>(e_idx, energies, particles, local_idx, particle_H, particle_G, dt);
    }

    static void accumulate(
        unsigned e_idx,
        const EnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat6r& particle_H,
        Vec6r& particle_G,
        Real dt
    )
    {
        accumulate<6>(e_idx, energies, particles, local_idx, particle_H, particle_G, dt);
    }
};

} // namespace Energy