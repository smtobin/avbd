#pragma once

#include "energy/TriangleRigidCollisionEnergyPool.hpp"
#include "energy/HardConstraintEnergySolver.hpp"

namespace Energy
{

struct TriangleRigidCollisionConstraintSolver
{
    /** Public typedefs */
    using PoolType = TriangleRigidCollisionEnergyPool;

    static Real evaluateConstraint(
        unsigned c_idx,
        const TriangleRigidCollisionEnergyPool& energies,
        ParticlePool& particles
    )
    {
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec4u& indices = energies.data[c_idx].particle_indices;

        const Vec3r& t1 = particles.positions[indices[0]];
        const Vec3r& t2 = particles.positions[indices[1]];
        const Vec3r& t3 = particles.positions[indices[2]];
        const Vec3r& barys = energies.data[c_idx].barys;

        const Vec3r cp_tri = barys[0]*t1 + barys[1]*t2 + barys[2]*t3;

        const Vec3r& rb_pos = particles.positions[indices[3]];
        const Quaternion& rb_rot = particles.rotation(indices[3]);
        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3r cp_rb = rb_pos + rb_rot * cp_rb_local;

        return n.dot(cp_rb - cp_tri);
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const TriangleRigidCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Real& C,
        Vec3r& C_grad,
        Mat3r& C_hess
    )
    {
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec4u& indices = energies.data[c_idx].particle_indices;

        if (local_idx < 3)
        {
            unsigned bary = energies.data[c_idx].barys[local_idx];
            C_grad = -bary * n;
            C_hess = Mat3r::Zero();
        }
        else
        {
            /** TODO: (07/21/26) gradient w.r.t oriented particle */
            C_grad = Vec3r::Zero();
            C_hess = Mat3r::Zero();
        }
    }
};

struct TriangleRigidCollisionEnergySolver
    : HardConstraintEnergySolver<TriangleRigidCollisionEnergyPool, TriangleRigidCollisionConstraintSolver>
{
    using HardConstraintEnergySolver::HardConstraintEnergySolver;
};

} // namespace Energy