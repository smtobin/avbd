#pragma once

#include "energy/GroundCollisionEnergyPool.hpp"
#include "energy/CollisionConstraintEnergySolver.hpp"

namespace Energy
{

struct GroundCollisionConstraintSolver
{
    /** Public typedefs */
    using PoolType = GroundCollisionEnergyPool;
    
    static void evaluateConstraint(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        Real& C_n, Real& C_t, Real& C_b
    )
    {
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        C_n = -particles.positions[p_idx][1];
        C_t = particles.positions[p_idx][0] - energies.data[c_idx].cp_x;
        C_b = particles.positions[p_idx][2] - energies.data[c_idx].cp_z;
    }

    template <int DOF>
    static void constraintGradientHessian(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned /* local_idx */,
        Real& C_n, Real& C_t, Real& C_b,
        Vec3r_or_Vec6r<DOF>& C_grad_n, Vec3r_or_Vec6r<DOF>& C_grad_t, Vec3r_or_Vec6r<DOF>& C_grad_b,
        Mat3r_or_Mat6r<DOF>& C_hess_n, Mat3r_or_Mat6r<DOF>& C_hess_t, Mat3r_or_Mat6r<DOF>& C_hess_b
    )
    {
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        C_n = -particles.positions[p_idx][1];
        C_t = particles.positions[p_idx][0] - energies.data[c_idx].cp_x;
        C_b = particles.positions[p_idx][2] - energies.data[c_idx].cp_z;

        // if (C > 0)
        //     particles.in_collision[p_idx] = true;


        if constexpr (DOF == 6)
        {
            C_grad_n = Vec6r(0, -1, 0, 0, 0, 0);
            C_grad_t = Vec6r(1, 0, 0, 0, 0, 0);
            C_grad_b = Vec6r(0, 0, 1, 0, 0, 0);
        }
        else
        {
            C_grad_n = Vec3r(0,-1,0);
            C_grad_t = Vec3r(1,0,0);
            C_grad_b = Vec3r(0,0,1);
        }
        


        C_hess_n = C_hess_t = C_hess_b = Mat3r_or_Mat6r<DOF>::Zero();
    }
};

// using GroundCollisionEnergySolver = 
//     HardConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>;

struct GroundCollisionEnergySolver
    : CollisionConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>
{
    using CollisionConstraintEnergySolver::CollisionConstraintEnergySolver;
    static constexpr bool SupportsPositional = true;
    static constexpr bool SupportsOriented = true;

    /** Updates the stiffness and Lagrange multiplier after a full time step.
     * Implements equation (19) from the AVBD paper.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     */
    static void updateAfterTimeStep(
        unsigned c_idx,
        GroundCollisionEnergyPool& energies,
        ParticlePool& particles
    )
    {
        // update the contact point
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        energies.data[c_idx].cp_x = particles.positions[p_idx][0];
        energies.data[c_idx].cp_z = particles.positions[p_idx][2];

        CollisionConstraintEnergySolver::updateAfterTimeStep(c_idx, energies, particles);
    }
};

} // namespace Energy