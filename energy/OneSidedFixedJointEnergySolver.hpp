#pragma once

#include "energy/OneSidedFixedJointEnergyPool.hpp"
#include "energy/HardConstraintEnergySolver.hpp"

#include "common/ParticlePool.hpp"
#include "common/Math.hpp"

namespace Energy
{

struct OneSidedFixedJointConstraintSolver
{
    using PoolType = OneSidedFixedJointEnergyPool;

    static Vec6r evaluateConstraint(
        unsigned c_idx,
        const OneSidedFixedJointEnergyPool& energies,
        ParticlePool& particles
    )
    {
        const OneSidedFixedJointEnergyInfo& info = energies.data[c_idx];
        const Vec3r& p = particles.positions[info.particle_indices[0]];
        const Quaternion& q = particles.rotation(info.particle_indices[0]);
        const Vec3r& jnt_p = p + q * info.position_offset;
        const Quaternion& jnt_q = q * info.rotation_offset;

        Vec6r C;
        C.head<3>() = jnt_p - info.ref_position;
        C.tail<3>() = Math::Minus_S3(jnt_q, info.ref_rotation);
        return C;

        std::cout << "C: " << C.transpose() << std::endl;
    }

    template <int>
    static void constraintGradientHessian(
        unsigned c_idx,
        const OneSidedFixedJointEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Vec6r& C,
        std::array<Vec6r, 6>& C_grad,
        std::array<Mat6r, 6>& C_hess
    )
    {
        const OneSidedFixedJointEnergyInfo& info = energies.data[c_idx];
        const Vec3r& p = particles.positions[info.particle_indices[0]];
        const Quaternion& q = particles.rotation(info.particle_indices[0]);
        const Vec3r& jnt_p = p + q * info.position_offset;
        const Quaternion& jnt_q = q * info.rotation_offset;

        Vec3r theta = Math::Minus_S3(jnt_q, info.ref_rotation);

        C.head<3>() = jnt_p - info.ref_position;
        C.tail<3>() = theta;

        std::cout << "C: " << C.transpose() << std::endl;

        for (int ci = 0; ci < 6; ci++)
        {
            C_grad[ci] = Vec6r::Zero();
            C_hess[ci] = Mat6r::Zero();
        }

        // gradient of positional DOF wrt position = Identity
        C_grad[0][0] = 1;
        C_grad[1][1] = 1;
        C_grad[2][2] = 1;


        // gradient of positional DOF wrt rotation
        Mat3r R = q.toRotationMatrix();
        Mat3r dCp_dR = -R * Math::Skew3(info.position_offset);
        C_grad[0].tail<3>() = dCp_dR.row(0);
        C_grad[1].tail<3>() = dCp_dR.row(1);
        C_grad[2].tail<3>() = dCp_dR.row(2);

        // gradient of rotational DOF wrt rotation
        Mat3r gam_inv = Math::ExpMap_InvRightJacobian(theta);
        Mat3r R_offset = info.rotation_offset.toRotationMatrix();
        Mat3r dCR_dR = gam_inv * R_offset.transpose();
        C_grad[3].tail<3>() = dCR_dR.row(0);
        C_grad[4].tail<3>() = dCR_dR.row(1);
        C_grad[5].tail<3>() = dCR_dR.row(2);


        /** TODO: (08/25/26) Hessian? Right now, just say it's 0. */

    }
};

struct OneSidedFixedJointEnergySolver
    : HardConstraintEnergySolver<OneSidedFixedJointEnergyPool, OneSidedFixedJointConstraintSolver>
{
    using Base = HardConstraintEnergySolver<OneSidedFixedJointEnergyPool, OneSidedFixedJointConstraintSolver>;
    using Base::Base;
    static constexpr bool SupportsPositional = false;
    static constexpr bool SupportsOriented = true;
};

    
} // namespace Energy