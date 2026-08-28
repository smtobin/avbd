#pragma once

#include "energy/collision/GroundCollisionEnergyPool.hpp"
#include "energy/collision/CollisionConstraintEnergySolver.hpp"

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

    /** Updates the contact point for this energy.
     * Since this energy is static (created upon simulation initialization), it is not updated by collision detection.
     * So, we must update the contact point here - this is called at the start of each time step.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     * @param particles : the particle pool
     * @param sdf_pool : the SDF primitive pool
     */
    static void updateContactPoints(
        unsigned c_idx,
        GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        Collision::SDFPrimitivePool& /* sdf_pool */
    )
    {
        // update the contact point
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        energies.data[c_idx].cp_x = particles.positions[p_idx][0];
        energies.data[c_idx].cp_z = particles.positions[p_idx][2];
    }
};

} // namespace Energy