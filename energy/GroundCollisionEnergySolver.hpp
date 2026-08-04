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
        C_t = particles.positions[p_idx][0];
        C_b = particles.positions[p_idx][2];
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Real& C_n, Real& C_t, Real& C_b,
        Vec3r& C_grad_n, Vec3r& C_grad_t, Vec3r& C_grad_b,
        Mat3r& C_hess_n, Mat3r& C_hess_t, Mat3r& C_hess_b
    )
    {
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        C_n = -particles.positions[p_idx][1];
        C_t = particles.positions[p_idx][0];
        C_b = particles.positions[p_idx][2];

        // if (C > 0)
        //     particles.in_collision[p_idx] = true;


        C_grad_n = Vec3r(0,-1,0);
        C_grad_t = Vec3r(1,0,0);
        C_grad_b = Vec3r(0,0,1);
        C_hess_n = C_hess_t = C_hess_b = Mat3r::Zero();
    }
};

// using GroundCollisionEnergySolver = 
//     HardConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>;

struct GroundCollisionEnergySolver
    : CollisionConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>
{
    using CollisionConstraintEnergySolver::CollisionConstraintEnergySolver;
};

} // namespace Energy