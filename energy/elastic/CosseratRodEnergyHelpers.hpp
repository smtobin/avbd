#pragma once

#include "common/common.hpp"
#include "common/Math.hpp"

namespace Energy
{

HD inline void cosseratRodEnergyAccumulate(
    const Vec3r& p1,
    const Vec3r& p2,
    const Quaternion& q1,
    const Quaternion& q2,
    Real rest_length,
    const Vec3r& precurvature,
    const Vec6r& stiffness,
    unsigned local_idx,
    Mat6r& H,
    Vec6r& G
)
{
    Real s_hat = 0.5;
    Real inv_length = 1.0/rest_length;

    Vec3r q_diff = Math::Minus_S3(q2, q1);
    Vec6r strain;

    // shear strain
    Quaternion q_mid = Math::Plus_S3(q1, s_hat * q_diff);
    strain.head<3>() = inv_length * (q_mid.conjugate() * (p2 - p1)) - Vec3r(0,0,1);

    // bending strain
    strain.tail<3>() = inv_length * q_diff - precurvature;

    Mat3r gam_inv = Math::ExpMap_InvRightJacobian(q_diff);

    // theta(s)
    Vec3r theta = s_hat * q_diff;

    // precompute useful quantities for gradients of shear strain
    Mat3r exp_theta = Math::Exp_so3(theta);
    Mat3r gam_theta = Math::ExpMap_RightJacobian(theta);
    Mat3r R = q_mid.toRotationMatrix();


    Vec3r dp_ds = inv_length * (p2 - p1);
    Mat3r RT_dp_ds_R = Math::Skew3(R.transpose() * dp_ds);

    Mat3r dtheta_dRi, dtheta_ds_dRi;
    Mat6r strain_grad;
    strain_grad.block<3,3>(3,0) = Mat3r::Zero();
    if (local_idx == 0)
    {
        dtheta_dRi = -s_hat * gam_inv.transpose();
        dtheta_ds_dRi = -inv_length * gam_inv.transpose();

        strain_grad.block<3,3>(0,0) = -inv_length * R.transpose();
        strain_grad.block<3,3>(0,3) =  RT_dp_ds_R * (exp_theta.transpose() + gam_theta * dtheta_dRi);
        strain_grad.block<3,3>(3,3) = dtheta_ds_dRi;
    }
    else
    {
        dtheta_dRi = s_hat * gam_inv;
        dtheta_ds_dRi = inv_length * gam_inv;

        strain_grad.block<3,3>(0,0) = inv_length * R.transpose();
        strain_grad.block<3,3>(0,3) = RT_dp_ds_R * gam_theta * dtheta_dRi;
        strain_grad.block<3,3>(3,3) = dtheta_ds_dRi;
    }

    // std::cout << "\nParticle indices: " << indices.transpose() << ", local_idx=" << local_idx << std::endl;
    // std::cout << "q1: \n" << q1 << "\nq2: \n" << q2 << std::endl;
    // std::cout << "q_mid: \n" << q_mid << "\nR mid: \n" << R << std::endl;
    // std::cout << "Strain: " << strain.transpose() << std::endl;
    // std::cout << "Strain gradient: \n" << strain_grad << std::endl;
    // std::cout << "Stiffness: " << stiffness.transpose() << std::endl;

    G = rest_length * strain_grad.transpose() * stiffness.cwiseProduct(strain);
    /** TODO: (08/23/26) Is Gauss-Newton approximation good enough? */
    Mat6r weighted_strain_grad = strain_grad.scaleRows(stiffness);
    Mat6r hess_an = rest_length * strain_grad.transpose() * weighted_strain_grad;
    // std::cout << "Hess diff: \n" << hess_fd - hess_an << std::endl;
    H = hess_an;
}

} // namespace Energy