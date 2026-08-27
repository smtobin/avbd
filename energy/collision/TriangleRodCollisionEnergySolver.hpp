#pragma once

#include "energy/collision/TriangleRodCollisionEnergyPool.hpp"
#include "energy/collision/CollisionConstraintEnergySolver.hpp"

namespace Energy
{

struct TriangleRodCollisionConstraintSolver
{
    /** Public typedefs */
    using PoolType = TriangleRodCollisionEnergyPool;

    static void evaluateConstraint(
        unsigned c_idx,
        const TriangleRodCollisionEnergyPool& energies,
        ParticlePool& particles,
        Real& C_n, Real& C_t, Real& C_b
    )
    {
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec3r& t = energies.data[c_idx].tangent;
        const Vec3r& b = energies.data[c_idx].binormal;
        const Vec5u& indices = energies.data[c_idx].particle_indices;

        const Vec3r& t1 = particles.positions[indices[0]];
        const Vec3r& t2 = particles.positions[indices[1]];
        const Vec3r& t3 = particles.positions[indices[2]];
        const Vec3r& barys = energies.data[c_idx].barys;

        const Vec3r cp_tri = barys[0]*t1 + barys[1]*t2 + barys[2]*t3;

        const Vec3r& sp1 = particles.positions[indices[3]];
        const Quaternion& sq1 = particles.rotation(indices[3]);
        const Vec3r& sp2 = particles.positions[indices[4]];
        const Quaternion& sq2 = particles.rotation(indices[4]);

        Real s = energies.data[c_idx].s;
        Quaternion q_mid = Math::Plus_S3(sq1, s*Math::Minus_S3(sq2, sq1));
        Vec3r p_mid = (1-s)*sp1 + s*sp2;

        const Vec3r& cp_rod_local = energies.data[c_idx].cp_rod_local;
        const Vec3r cp_seg = p_mid + q_mid*cp_rod_local;

        Vec3r diff = cp_seg - cp_tri;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);
    }

    template <int DOF>
    static void constraintGradientHessian(
        unsigned c_idx,
        const TriangleRodCollisionEnergyPool& energies,
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
        const Vec5u& indices = energies.data[c_idx].particle_indices;

        const Vec3r& t1 = particles.positions[indices[0]];
        const Vec3r& t2 = particles.positions[indices[1]];
        const Vec3r& t3 = particles.positions[indices[2]];
        const Vec3r& barys = energies.data[c_idx].barys;

        const Vec3r cp_tri = barys[0]*t1 + barys[1]*t2 + barys[2]*t3;

        const Vec3r& sp1 = particles.positions[indices[3]];
        const Quaternion& sq1 = particles.rotation(indices[3]);
        const Vec3r& sp2 = particles.positions[indices[4]];
        const Quaternion& sq2 = particles.rotation(indices[4]);

        Real s = energies.data[c_idx].s;
        Vec3r q_diff = Math::Minus_S3(sq2, sq1);
        Quaternion q_mid = Math::Plus_S3(sq1, s*q_diff);
        Vec3r p_mid = (1-s)*sp1 + s*sp2;

        const Vec3r& cp_rod_local = energies.data[c_idx].cp_rod_local;
        const Vec3r cp_seg = p_mid + q_mid*cp_rod_local;

        Vec3r diff = cp_seg - cp_tri;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);

        if constexpr (DOF == 6)
        {
            /** TODO: (07/21/26) gradient w.r.t oriented particle */
            Mat3r R = q_mid.toRotationMatrix();
            Mat3r skew_cp = Math::Skew3(cp_rod_local);
            Mat3r R_rloc = -R * skew_cp;

            Mat3r gam_inv = Math::ExpMap_InvRightJacobian(q_diff);
            Mat3r gam = Math::ExpMap_RightJacobian(s*q_diff);
            

            Mat3r dRs_dR;
            Real scaling;
            if (local_idx == 3)
            {
                Mat3r exp = Math::Exp_so3(s*q_diff);
                dRs_dR = exp.transpose() - s * gam * gam_inv.transpose();

                scaling = 1-s;
            }
            else if (local_idx == 4)
            {
                dRs_dR = s * gam * gam_inv;

                scaling = s;
            }
            else
                throw std::runtime_error("TriangleRigidCollisionEnergySolver::constraintGradientHessian local_idx out of bounds!");
            
            Mat3r R_rloc_dRs_dR = R_rloc * dRs_dR;
            
            // gradients
            C_grad_n.template block<3,1>(0,0) = scaling*n;
            C_grad_n.template block<3,1>(3,0) = -n.transpose() * R_rloc_dRs_dR;
            C_grad_t.template block<3,1>(0,0) = scaling*t;
            C_grad_t.template block<3,1>(3,0) = -t.transpose() * R_rloc_dRs_dR;
            C_grad_b.template block<3,1>(0,0) = scaling*b;
            C_grad_b.template block<3,1>(3,0) = -b.transpose() * R_rloc_dRs_dR;
            C_grad_t = C_grad_b = Vec6r::Zero();

            // Hessians set to 0 (for now just use Gauss-Newton)
            C_hess_n = C_hess_t = C_hess_b = Mat6r::Zero();
            

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

struct TriangleRodCollisionEnergySolver
    : CollisionConstraintEnergySolver<TriangleRodCollisionEnergyPool, TriangleRodCollisionConstraintSolver>
{
    using CollisionConstraintEnergySolver::CollisionConstraintEnergySolver;
    static constexpr bool SupportsPositional = true;
    static constexpr bool SupportsOriented = true;
};

} // namespace Energy