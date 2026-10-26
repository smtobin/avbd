#pragma once

#include "common/common.hpp"

namespace Energy
{

/** Helper function to compute the deformation gradient F given 4 vertices and the rest-state matrix Q.
 * 
 * F = [v1 - v4   v2 - v4   v3 - v4] * Q
 */
HD inline void computeDeformationGradient(
    const Vec3r& v1,
    const Vec3r& v2,
    const Vec3r& v3, 
    const Vec3r& v4,
    const Mat3r& Q,
    Mat3r& F
)
{
    F.setCol(0, v1 - v4);
    F.setCol(1, v2 - v4);
    F.setCol(2, v3 - v4);

    F = F*Q;
}

HD inline void neoHookeanEnergyAccumulate(
    const Vec3r& p0,
    const Vec3r& p1,
    const Vec3r& p2,
    const Vec3r& p3,
    Real V,
    const Mat3r& Q,
    Real lambda,
    Real mu,
    Real kd,
    const Mat3r& E_prev,
    Real dt,
    unsigned local_idx,
    Mat3r& H,
    Vec3r& G
    
)
{
    // compute F for this timestep and the previous timestep
    Mat3r F;
    computeF(p0, p1, p2, p3, Q, F);

    // hydrostatic gradient
    Vec3r c0 = F.col(1).cross(F.col(2));
    Vec3r c1 = F.col(2).cross(F.col(0));
    Vec3r c2 = F.col(0).cross(F.col(1));

    const Real detF = F.col(0).dot(c0);

    // compute the Hessian and gradient
    Vec3r qi;
    Real sign;
    if (local_idx < 3)
    {
        qi = Q.row(local_idx);
        sign = 1.0;
    }
    else
    {
        qi = Q.row(0) + Q.row(1) + Q.row(2);
        sign = -1.0;
    }

    const Vec3r Fqi = F*qi;
    const Vec3r Fcross_qi = c0*qi[0] + c1*qi[1] + c2*qi[2];
    
    const Real hyd_mult = V*lambda * (detF - mu/lambda - 1);
    const Real dev_mult = V*mu;

    // compute rate of change of Green strain
    Mat3r FtF = F.transpose() * F;
    
    const Real qi_norm2 = qi.squaredNorm();

    // gradient
    Real damp_mult_factor = (2/dt * V * kd);
    G = sign * hyd_mult * Fcross_qi;
    G += sign * dev_mult * Fqi;

    Mat3r E_dot = 0.5 * FtF;
    E_dot(0,0) -= 0.5; E_dot(1,1) -= 0.5; E_dot(2,2) -= 0.5;
    E_dot -= E_prev;
    E_dot *= 1/dt;

    G += (sign * damp_mult_factor) * (F * (E_dot * qi));

    // Hessian
    // hydrostatic
    H = (V * lambda) * outerProduct(Fcross_qi, Fcross_qi);
    // deviatoric
    Real t1 = V * mu * qi_norm2;
    H(0,0) += t1; H(1,1) += t1; H(2,2) += t1;

    // damping
    Real qi_Edot_qi = qi.dot(E_dot * qi);
    Real t2 = qi_Edot_qi * damp_mult_factor;
    H(0,0) += t2; H(1,1) += t2; H(2,2) += t2;

    H += (0.5 * damp_mult_factor/dt) * (F*F.transpose() * qi_norm2 + outerProduct(Fqi, Fqi));
}

}