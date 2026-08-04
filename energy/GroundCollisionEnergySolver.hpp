#pragma once

#include "energy/GroundCollisionEnergyPool.hpp"
#include "energy/CollisionConstraintEnergySolver.hpp"

namespace Energy
{

struct GroundCollisionConstraintSolver
{
    /** Public typedefs */
    using PoolType = GroundCollisionEnergyPool;
    
    static Vec3r evaluateConstraint(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles
    )
    {
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        Real C_n = -particles.positions[p_idx][1];
        Real C_t = particles.positions[p_idx][0];
        Real C_b = particles.positions[p_idx][2];

        return Vec3r(C_t, C_b, C_n);
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Vec3r& C,
        Vec3r& C_grad_n, Vec3r& C_grad_t, Vec3r& C_grad_b,
        Mat3r& C_hess_n, Mat3r& C_hess_t, Mat3r& C_hess_b
    )
    {
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        C[2] = -particles.positions[p_idx][1];
        C[0] = particles.positions[p_idx][0];
        C[1] = particles.positions[p_idx][2];

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