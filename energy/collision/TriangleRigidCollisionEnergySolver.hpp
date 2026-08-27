#pragma once

#include "energy/collision/TriangleRigidCollisionEnergyPool.hpp"
#include "energy/collision/CollisionConstraintEnergySolver.hpp"

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

    template <int DOF>
    static void constraintGradientHessian(
        unsigned c_idx,
        const TriangleRigidCollisionEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Real& C_n, Real& C_t, Real& C_b,
        Vec3r_or_Vec6r<DOF>& C_grad_n, Vec3r_or_Vec6r<DOF>& C_grad_t, Vec3r_or_Vec6r<DOF>& C_grad_b,
        Mat3r_or_Mat6r<DOF>& C_hess_n, Mat3r_or_Mat6r<DOF>& C_hess_t, Mat3r_or_Mat6r<DOF>& C_hess_b
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

        if constexpr (DOF == 6)
        {
            /** TODO: (07/21/26) gradient w.r.t oriented particle */
            Mat3r R = rb_rot.toRotationMatrix();
            Mat3r skew_cp = Math::Skew3(cp_rb_local);
            Mat3r R_rloc = R * skew_cp;
            // gradients
            C_grad_n.template block<3,1>(0,0) = n;
            C_grad_n.template block<3,1>(3,0) = -n.transpose() * R_rloc;
            C_grad_t.template block<3,1>(0,0) = t;
            C_grad_t.template block<3,1>(3,0) = -t.transpose() * R_rloc;
            C_grad_b.template block<3,1>(0,0) = b;
            C_grad_b.template block<3,1>(3,0) = -b.transpose() * R_rloc;
            C_grad_t = C_grad_b = Vec6r::Zero();

            // Hessians
            // grad = (Skew(cp_rb_local) * R^T * n)^T
            // ==> (Skew(cp_rb_local) * -R^T * skew(n) * -R)^T = R^T * skew(n) * R * skew(cp_rb_local)
            C_hess_n = C_hess_t = C_hess_b = Mat6r::Zero();
            C_hess_n.template block<3,3>(3,3) = R.transpose() * Math::Skew3(n) * R * skew_cp;
            C_hess_n.template block<3,3>(3,3) = 0.5*(C_hess_n.template block<3,3>(3,3) + C_hess_n.template block<3,3>(3,3).transpose());
            C_hess_t.template block<3,3>(3,3) = R.transpose() * Math::Skew3(t) * R * skew_cp;
            C_hess_t.template block<3,3>(3,3) = 0.5*(C_hess_t.template block<3,3>(3,3) + C_hess_t.template block<3,3>(3,3).transpose());
            C_hess_b.template block<3,3>(3,3) = R.transpose() * Math::Skew3(b) * R * skew_cp;
            C_hess_b.template block<3,3>(3,3) = 0.5*(C_hess_b.template block<3,3>(3,3) + C_hess_b.template block<3,3>(3,3).transpose());

        }
        else
        {
            if (local_idx >= 3)
                throw std::runtime_error("TriangleRigidCollisionEnergySolver::constraintGradientHessian local_idx out of bounds!");

            Real bary = energies.data[c_idx].barys[local_idx];
            C_grad_n = -bary * n;
            C_grad_t = -bary * t;
            C_grad_b = -bary * b;
            C_hess_n = Mat3r::Zero();
            C_hess_t = Mat3r::Zero();
            C_hess_b = Mat3r::Zero();
        }
    }
};

struct TriangleRigidCollisionEnergySolver
    : CollisionConstraintEnergySolver<TriangleRigidCollisionEnergyPool, TriangleRigidCollisionConstraintSolver>
{
    using CollisionConstraintEnergySolver::CollisionConstraintEnergySolver;
    static constexpr bool SupportsPositional = true;
    static constexpr bool SupportsOriented = true;
};

} // namespace Energy