#pragma once

#include "energy/TriangleRigidCollisionEnergyPool.hpp"
#include "energy/CollisionConstraintEnergySolver.hpp"

namespace Energy
{

struct TriangleRigidCollisionConstraintSolver
{
    /** Public typedefs */
    using PoolType = TriangleRigidCollisionEnergyPool;

    static void evaluateConstraint(
        unsigned c_idx,
        const TriangleRigidCollisionEnergyPool& energies,
        ParticlePool& particles,
        Real& C_n, Real& C_t, Real& C_b
    )
    {
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec3r& t = energies.data[c_idx].tangent;
        const Vec3r& b = energies.data[c_idx].binormal;
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

        Vec3r diff = cp_rb - cp_tri;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);

        // if (C > 0)
        // {
        //     particles.in_collision[indices[0]] = true;
        //     particles.in_collision[indices[1]] = true;
        //     particles.in_collision[indices[2]] = true;
        //     particles.in_collision[indices[3]] = true;
        // }
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const TriangleRigidCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Real& C_n, Real& C_t, Real& C_b,
        Vec3r& C_grad_n, Vec3r& C_grad_t, Vec3r& C_grad_b,
        Mat3r& C_hess_n, Mat3r& C_hess_t, Mat3r& C_hess_b
    )
    {
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec3r& t = energies.data[c_idx].tangent;
        const Vec3r& b = energies.data[c_idx].binormal;
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

        Vec3r diff = cp_rb - cp_tri;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);

        // if (C > 0)
        // {
        //     particles.in_collision[indices[0]] = true;
        //     particles.in_collision[indices[1]] = true;
        //     particles.in_collision[indices[2]] = true;
        //     particles.in_collision[indices[3]] = true;
        // }

        if (local_idx < 3)
        {
            Real bary = energies.data[c_idx].barys[local_idx];
            C_grad_n = -bary * n;
            C_grad_t = -bary * t;
            C_grad_b = -bary * b;
            C_hess_n = Mat3r::Zero();
            C_hess_t = Mat3r::Zero();
            C_hess_b = Mat3r::Zero();
        }
        else
        {
            /** TODO: (07/21/26) gradient w.r.t oriented particle */
            C_grad_n = C_grad_t = C_grad_b = Vec3r::Zero();
            C_hess_n = C_hess_t = C_hess_b = Mat3r::Zero();
        }
    }
};

struct TriangleRigidCollisionEnergySolver
    : CollisionConstraintEnergySolver<TriangleRigidCollisionEnergyPool, TriangleRigidCollisionConstraintSolver>
{
    using CollisionConstraintEnergySolver::CollisionConstraintEnergySolver;
};

} // namespace Energy