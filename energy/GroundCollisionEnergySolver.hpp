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
        unsigned particle_idx = energies.particles[c_idx];
        return -particles.positions[particle_idx][1];
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const GroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        Real& C,
        Vec3r grads[1],
        Mat3r hessians[1]
    )
    {
        unsigned p_idx = energies.particles[c_idx];
        C = -particles.positions[p_idx][1];

        grads[0] = Vec3r(0,-1,0);
        hessians[0] = Mat3r::Zero();
    }
};

using GroundCollisionEnergySolver = 
    HardConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>;

} // namespace Energy