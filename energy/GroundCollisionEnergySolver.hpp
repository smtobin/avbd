#pragma once

#include "energy/GroundCollisionEnergyPool.hpp"
#include "energy/HardConstraintEnergySolver.hpp"

namespace Energy
{

struct GroundCollisionConstraintSolver
{
    static Real evaluateConstraint(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles
    )
    {
        unsigned particle_idx = energies.data[c_idx].particle_indices[0];
        return -particles.positions[particle_idx][1];
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Real& C,
        Vec3r& C_grad,
        Mat3r& C_hess
    )
    {
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        C = -particles.positions[p_idx][1];

        C_grad = Vec3r(0,-1,0);
        C_hess = Mat3r::Zero();
    }
};

// using GroundCollisionEnergySolver = 
//     HardConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>;

struct GroundCollisionEnergySolver
    : HardConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>
{
    using HardConstraintEnergySolver::HardConstraintEnergySolver;
};

} // namespace Energy