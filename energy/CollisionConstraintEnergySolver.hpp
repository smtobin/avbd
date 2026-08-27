#pragma once

#include "CollisionConstraintEnergyPool.hpp"
#include "common/Math.hpp"

namespace Energy
{

template <typename CollisionEnergyPool, typename CollisionConstraintSolver>
struct CollisionConstraintEnergySolver
{
    using PoolType = CollisionEnergyPool;
    
    /** Updates the collision normal, which causes an update to the collision tangent and binormal.
     * The Lagrange multipliers associated with the tangent and binormal are updated so that the total force remains the same.
     */
    static void updateCollisionFrame(
        unsigned c_idx,
        CollisionEnergyPool& energies,
        const Vec3r& new_normal
    )
    {
        // extract old n, t, b
        Vec3r& t = energies.data[c_idx].tangent;
        Vec3r& b = energies.data[c_idx].binormal;
        Vec3r& n = energies.data[c_idx].normal;

        Real& lambda_t = energies.data[c_idx].lambda_t;
        Real& lambda_b = energies.data[c_idx].lambda_b;

        // compute old frictional force
        Vec3r f_tb_old = lambda_t * t + lambda_b * b;

        // project last frame's tangent vector onto the plane orthogonal to the new normal vector
        Vec3r t_proj = t - (t.dot(new_normal))*new_normal;
        // renormalize - this is the new tangent vector given the new normal
        // designed to be coherent from previous frames, so that the collision frame does not change too drastically between frames
        Real t_proj_mag = t_proj.norm();
        if (t_proj_mag > 1e-8)
        {
            t = t_proj / t_proj_mag;
            b = new_normal.cross(t);
        }
        else
        {
            Math::completeOrthonormalBasisGivenNormal(new_normal, t, b);
        }
        n = new_normal;

        // update tangent and binormal lambda so that friction force stays the same
        lambda_t = f_tb_old.dot(t);
        lambda_b = f_tb_old.dot(b);
    }

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
        Real C_n, C_t, C_b;
        CollisionConstraintSolver::evaluateConstraint(c_idx, energies, particles, C_n, C_t, C_b);

        // subtract off previous constraint violation
        Real C_corr_n = C_n - CONSTRAINT_ALPHA * energies.data[c_idx].C_n_prev;
        
        // extract the current stiffnesses and Lagrange multipliers
        Real& k_n = energies.data[c_idx].k_n;
        Real& k_t = energies.data[c_idx].k_t;
        Real& k_b = energies.data[c_idx].k_b;
        Real& lambda_n = energies.data[c_idx].lambda_n;
        Real& lambda_t = energies.data[c_idx].lambda_t;
        Real& lambda_b = energies.data[c_idx].lambda_b;
        bool& use_static = energies.data[c_idx].use_static;

        Real mu_s = energies.data[c_idx].mu_s;
        Real mu_k = energies.data[c_idx].mu_k;

        Real lambda_n_plus = k_n * C_corr_n + lambda_n;
        Real lambda_t_plus = k_t * C_t + lambda_t;
        Real lambda_b_plus = k_b * C_b + lambda_b;

        // update normal stiffness - equation (12)
        // for collision constraints, lambda min = 0
        if (lambda_n_plus > 0)
        {
            // additive stiffness update
            k_n += STIFFNESS_BETA * std::abs(C_corr_n);
        }

        // update normal lambda - equation (11)
        lambda_n = std::max(Real(0), lambda_n_plus);

        // update tangent and binormal stiffness - equation (12)
        Vec2r lambda_tb_plus(lambda_t_plus, lambda_b_plus);
        Real lambda_tb_plus_mag = lambda_tb_plus.norm();
        Real mu = use_static ? mu_s : mu_k;
        Real lambda_tb_max = mu*lambda_n;
        if (lambda_tb_plus_mag < lambda_tb_max)
        {
            k_t += STIFFNESS_BETA * std::abs(C_t);
            k_b += STIFFNESS_BETA * std::abs(C_b);

            lambda_t = lambda_t_plus;
            lambda_b = lambda_b_plus;

            energies.data[c_idx].unclamped_lambda_tb = true;
        }
        else
        {
            // if we must clamp the force and static friction was used, we have overcome the static friction boundary
            // so set use_static to false, and immediatedly update the max lambda_tb allowed to use kinetic friction coefficient
            if (use_static)
            {
                use_static = false;
                lambda_tb_max = mu_k * lambda_n_plus;
            }

            Vec2r lambda_tb_clamped = lambda_tb_plus / (lambda_tb_plus_mag + 1e-8) * lambda_tb_max;
            lambda_t = lambda_tb_clamped[0];
            lambda_b = lambda_tb_clamped[1];

            
            energies.data[c_idx].unclamped_lambda_tb = false;
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
        ParticlePool& particles
    )
    {
        // store the current constraint violation
        if (energies.data[c_idx].lambda_n > 0)
        {
            Real C_n, C_t, C_b;
            CollisionConstraintSolver::evaluateConstraint(c_idx, energies, particles, C_n, C_t, C_b);
            energies.data[c_idx].C_n_prev = C_n;
        }

        // if, at the final iteration of the time step, the frictional Lagrange multiplier magnitude was unclamped
        // (i.e. friction erased all tangential movement for this particle)
        // then if we were using kinetic friction, switch to static
        if (energies.data[c_idx].unclamped_lambda_tb && !energies.data[c_idx].use_static)
            energies.data[c_idx].use_static = true;

        energies.data[c_idx].lambda_n *= CONSTRAINT_ALPHA * STIFFNESS_GAMMA;
        energies.data[c_idx].lambda_t *= STIFFNESS_GAMMA;
        energies.data[c_idx].lambda_b *= STIFFNESS_GAMMA;
        energies.data[c_idx].k_n = std::max(STIFFNESS_GAMMA * energies.data[c_idx].k_n, energies.k_start);
        energies.data[c_idx].k_t = std::max(STIFFNESS_GAMMA * energies.data[c_idx].k_t, energies.k_start);
        energies.data[c_idx].k_b = std::max(STIFFNESS_GAMMA * energies.data[c_idx].k_b, energies.k_start);
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
        const CollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat3r_or_Mat6r<DOF>& particle_H,
        Vec3r_or_Vec6r<DOF>& particle_G,
        Real dt
    )
    {
        // evaluate the constraint, gradients, and Hessians for the constraint for each particle involved
        Real C_n, C_t, C_b;
        Vec3r_or_Vec6r<DOF> C_grad_n, C_grad_t, C_grad_b;
        Mat3r_or_Mat6r<DOF> C_hess_n, C_hess_t, C_hess_b;
        CollisionConstraintSolver::template constraintGradientHessian<DOF>(
            e_idx, energies, particles, local_idx, // inputs
            C_n, C_t, C_b, 
            C_grad_n, C_grad_t, C_grad_b,
            C_hess_n, C_hess_t, C_hess_b  // outputs
        );
        // subtract off previous constraint violation
        Real C_corr_n = C_n - CONSTRAINT_ALPHA * energies.data[e_idx].C_n_prev;

        // extract the current stiffnesses and Lagrange multipliers
        Real k_n = energies.data[e_idx].k_n;
        Real k_t = energies.data[e_idx].k_t;
        Real k_b = energies.data[e_idx].k_b;
        Real lambda_n = energies.data[e_idx].lambda_n;
        Real lambda_t = energies.data[e_idx].lambda_t;
        Real lambda_b = energies.data[e_idx].lambda_b;

        bool use_static = energies.data[e_idx].use_static;
        Real mu_s = energies.data[e_idx].mu_s;
        Real mu_k = energies.data[e_idx].mu_k;


        // Lagrange multiplier for normal
        Real lambda_n_plus_unclamped = k_n * C_corr_n + lambda_n;
        Real lambda_n_plus = std::max(Real(0), lambda_n_plus_unclamped);

        // Lagrange multipliers for tangent and binormal
        Vec2r lambda_tb_plus(k_t*C_t + lambda_t, k_b*C_b + lambda_b);
        Vec2r lambda_tb_plus_unclamped = lambda_tb_plus;
        // clamp magnitude
        Real lambda_tb_plus_mag = lambda_tb_plus.norm();
        Real mu = use_static ? mu_s : mu_k;
        Real lambda_tb_max = mu*lambda_n_plus;
        if (lambda_tb_plus_mag > lambda_tb_max)
        {
            lambda_tb_plus = lambda_tb_plus / (lambda_tb_plus.norm() + 1e-8) * lambda_tb_max;
        }
        
        
        /** TODO: (08/04/26) After we've clamped the lambdas, these will always be false...? */
        // stiffness rescaling for normal - equation (14)
        Real k_scaled_n = k_n;
        if (lambda_n_plus_unclamped < 0 && std::abs(C_n) > 1e-12)
            k_scaled_n =  std::abs(-lambda_n / C_n);
        
        // stiffness rescling for tangent and binormal - equation (14)
        Vec2r k_scaled_tb(k_t, k_b);
        Vec2r C_tb(C_t, C_b);
        Vec2r lambda_tb(lambda_t, lambda_b);
        if (lambda_tb_plus_mag > lambda_tb_max && C_tb.norm() > 1e-12)
        {
            k_scaled_tb = k_scaled_tb / k_scaled_tb.norm() * std::abs(lambda_tb_max - lambda_tb.norm()) / C_tb.norm();
        }
        // gradient
        Vec3r_or_Vec6r<DOF> grad_n = lambda_n_plus * C_grad_n;
        Vec3r_or_Vec6r<DOF> grad_t = lambda_tb_plus[0] * C_grad_t;
        Vec3r_or_Vec6r<DOF> grad_b = lambda_tb_plus[1] * C_grad_b;

            // Hessian
            /** TODO: do diagonalization of Hessian component
             * 
             * 
             * 
             * 
             */
            

        Mat3r_or_Mat6r<DOF> hess_n = (k_scaled_n * C_corr_n + lambda_n) * C_hess_n +
            k_scaled_n * C_grad_n* C_grad_n.transpose();
        Mat3r_or_Mat6r<DOF> hess_t = (k_scaled_tb[0] * C_t + lambda_t) * C_hess_t +
            k_scaled_tb[0] * C_grad_t * C_grad_t.transpose();
        Mat3r_or_Mat6r<DOF> hess_b = (k_scaled_tb[1] * C_b + lambda_b) * C_hess_b +
            k_scaled_tb[1] * C_grad_b * C_grad_b.transpose();

        particle_H += hess_n + hess_t + hess_b;
        particle_G += grad_n + grad_t + grad_b;
    }

    static void accumulate(
        unsigned e_idx,
        const CollisionEnergyPool& energies,
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