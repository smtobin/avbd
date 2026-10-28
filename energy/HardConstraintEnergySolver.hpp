#pragma once

#include "common/common.hpp"
#include "energy/HardConstraintEnergyPool.hpp"
#include "common/ParticlePool.hpp"


namespace Energy
{

template <typename EnergyPool, typename ConstraintSolver>
struct HardConstraintEnergySolver
{
    using Vec = VecNr<EnergyPool::Dim>;

    /** Energy function */
    static Real energy(
        unsigned e_idx,
        const typename EnergyPool::ViewType& energies,
        ParticlePoolView& particles,
        Real /* dt */
    )
    {
        // evaluate the constraint
        Vec C = ConstraintSolver::evaluateConstraint(e_idx, energies, particles);
        const Vec& k = energies.data[e_idx].k;
        const Vec& lambda = energies.data[e_idx].lambda;

        Real E = 0;
        for (unsigned ci = 0; ci < EnergyPool::Dim; ci++)
        {
            E += 0.5 * k[ci] * C[ci] * C[ci] + lambda[ci] * C[ci];
        }
        return E;
    }

    /** Updates the stiffness and Lagrange multiplier after the iteration
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     * @param particles : the simulation particle memory pool
     */
    static void updateAfterIteration(
        unsigned c_idx,
        typename EnergyPool::ViewType& energies,
        ParticlePoolView& particles
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
        typename EnergyPool::ViewType& energies,
        ParticlePoolView& particles
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
     * @param H : the current accumulated Hessian - added to by this function
     * @param G : the current accumulated gradient - added to by this function
     * @param dt : the time step
     */
    template <int DOF>
    static void accumulate(
        unsigned e_idx,
        const typename EnergyPool::ViewType& energies,
        ParticlePoolView& particles,
        unsigned local_idx,
        MatNr<DOF>& H,
        VecNr<DOF>& G,
        Real /* dt */
    )
    {
        std::cout << "HardconstraintEnergy accumulate! DOF=" << DOF << std::endl;
        // evaluate the constraint, gradients, and Hessians for the constraint for each particle involved
        std::array<VecNr<DOF>, EnergyPool::Dim> C_grad;
        std::array<MatNr<DOF>, EnergyPool::Dim> C_hess;       // constraint Hessian is 3rd order tensor
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
            G += C_grad[ci] * lambda_p[ci];

                // Hessian
                /** TODO: do diagonalization of Hessian component
                 * 
                 * 
                 * 
                 * 
                 */
                

            H += (k_scaled[ci] * C[ci] + lambda[ci]) * C_hess[ci] +
                k_scaled[ci] * outerProduct(C_grad[ci], C_grad[ci]);
        }
    }

    static void accumulate(
        unsigned e_idx,
        const typename EnergyPool::ViewType& energies,
        ParticlePoolView& particles,
        unsigned local_idx,
        Mat3r& H,
        Vec3r& G,
        Real dt
    )
    {
        accumulate<3>(e_idx, energies, particles, local_idx, H, G, dt);
    }

    static void accumulate(
        unsigned e_idx,
        const typename EnergyPool::ViewType& energies,
        ParticlePoolView& particles,
        unsigned local_idx,
        Mat6r& H,
        Vec6r& G,
        Real dt
    )
    {
        accumulate<6>(e_idx, energies, particles, local_idx, H, G, dt);
    }
};

} // namespace Energy