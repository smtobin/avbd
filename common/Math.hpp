#pragma once

#include "common/common.hpp"
#include "common/Quaternion.hpp"
#include <iostream>


class Math
{

public:

static Mat3r Skew3(const Vec3r& vec)
{
    Mat3r mat;
    mat << 0,       -vec(2),    vec(1),
           vec(2),  0,          -vec(0),
           -vec(1), vec(0),     0;
    return mat;
}

static Vec3r Vee3(const Mat3r& mat)
{
    return Vec3r(mat(2,1), mat(0,2), mat(1,0));
}

// Computes right Jacobian of SO(3) exponential map
static Mat3r ExpMap_RightJacobian(const Vec3r& theta)
{
    Real theta_norm = theta.norm();
    const Mat3r skew = Skew3(theta);
    if (theta_norm < Real(1e-8))
    {
        return Mat3r::Identity() - 0.5 * skew;
    }
    
    return Mat3r::Identity() - (1 - std::cos(theta_norm)) / (theta_norm * theta_norm) * skew + (theta_norm - std::sin(theta_norm)) / (theta_norm * theta_norm * theta_norm) * skew * skew;
}

// Computes derivative of right Jacobian of SO(3) exponential map (at theta), multiplied by some perturbation xi
static Mat3r DExpMap_RightJacobian_Contract_k(const Vec3r& theta, const Vec3r& xi)
{
    Real theta_norm = theta.norm();

    Mat3r skew_theta = Skew3(theta);
    Mat3r skew_xi = Skew3(xi);

    if (theta_norm < Real(1e-8))
    {
        return  1.0/12.0 * theta.dot(xi) * skew_theta - 0.5 * skew_xi;
    }

    Real theta_norm2 = theta_norm*theta_norm;
    Real theta_norm3 = theta_norm2*theta_norm;
    Real theta_norm4 = theta_norm3*theta_norm;

    Real a = (1 - std::cos(theta_norm)) / theta_norm2;
    Real b = (theta_norm - std::sin(theta_norm)) / theta_norm3;
    Real da_dtheta_norm = -2*(1 - std::cos(theta_norm)) / theta_norm3 + std::sin(theta_norm) / theta_norm2;
    Real db_dtheta_norm = (1 - std::cos(theta_norm)) / theta_norm3 - 3*(theta_norm - std::sin(theta_norm)) / theta_norm4;

    

    Mat3r term1 = -da_dtheta_norm * theta.dot(xi) / theta_norm * skew_theta;
    Mat3r term2 = -a * skew_xi;
    Mat3r term3 = db_dtheta_norm * theta.dot(xi) / theta_norm * skew_theta * skew_theta;
    Mat3r term4 = b * (skew_xi * skew_theta + skew_theta * skew_xi);

    return term1 + term2 + term3 + term4;
    
}

static Mat3r DExpMap_RightJacobian_Contract_j(const Vec3r& theta, const Vec3r& xi)
{
    Real theta_norm = theta.norm();

    if (theta_norm < Real(1e-2))
    {
        return DExpMap_RightJacobian_Contract_j_approx(theta, xi);
    }

    Mat3r skew_theta = Skew3(theta);

    Mat3r dskew_dtheta_i[3] = {Mat3r::Zero(), Mat3r::Zero(), Mat3r::Zero()};
    dskew_dtheta_i[0](1, 2) = -1;
    dskew_dtheta_i[0](2, 1) = 1;

    dskew_dtheta_i[1](0, 2) = 1;
    dskew_dtheta_i[1](2, 0) = -1;

    dskew_dtheta_i[2](0, 1) = -1;
    dskew_dtheta_i[2](1, 0) = 1;

    Real theta_norm2 = theta_norm*theta_norm;
    Real theta_norm3 = theta_norm2*theta_norm;
    Real theta_norm4 = theta_norm3*theta_norm;

    Real a = (1 - std::cos(theta_norm)) / theta_norm2;
    Real b = (theta_norm - std::sin(theta_norm)) / theta_norm3;
    Real da_dtheta_norm = -2*(1 - std::cos(theta_norm)) / theta_norm3 + std::sin(theta_norm) / theta_norm2;
    Real db_dtheta_norm = (1 - std::cos(theta_norm)) / theta_norm3 - 3*(theta_norm - std::sin(theta_norm)) / theta_norm4;


    Mat3r term1, term2, term3, term4;
    for (int k = 0; k < 3; k++)
    {
        term1.col(k) = -da_dtheta_norm * theta[k] / theta_norm * skew_theta * xi;
        term2.col(k) = -a * dskew_dtheta_i[k] * xi;
        term3.col(k) = db_dtheta_norm * theta[k] / theta_norm * skew_theta * skew_theta * xi;
        term4.col(k) = b * (dskew_dtheta_i[k] * skew_theta + skew_theta * dskew_dtheta_i[k]) * xi;
    }

    return term1 + term2 + term3 + term4;

}

static Mat3r DExpMap_RightJacobian_Contract_j_approx(const Vec3r& theta, const Vec3r& xi)
{
    Vec3r skew_xi = theta.cross(xi);

    Vec3r dskew_dtheta_xi_i[] = { Vec3r(0,-xi[2],xi[1]), Vec3r(xi[2],0,-xi[0]), Vec3r(-xi[1],xi[0],0) };
    Vec3r dskew_dtheta_skew_xi_i[] = { Vec3r(0,-skew_xi[2],skew_xi[1]), Vec3r(skew_xi[2],0,-skew_xi[0]), Vec3r(-skew_xi[1],skew_xi[0],0) };

    Mat3r term1, term2, term3;
    for (int k = 0; k < 3; k++)
    {
        term1.col(k) = 1.0/12.0 * theta[k] * skew_xi;
        term2.col(k) = -1.0/2.0 * dskew_dtheta_xi_i[k];
        term3.col(k) = 1.0/6.0 * (dskew_dtheta_skew_xi_i[k] + theta.cross(dskew_dtheta_xi_i[k]));
    }

    return term1 + term2 + term3;
}

// Computes the inverse of the right Jacobian of SO(3) exponential map
static Mat3r ExpMap_InvRightJacobian(const Vec3r& theta)
{
    Real theta_norm = theta.norm();
    const Mat3r skew = Skew3(theta);

    if (theta_norm < Real(1e-8))
    {
        return Mat3r::Identity() + 0.5 * skew;
    }

    return Mat3r::Identity() + 0.5*skew + (1/(theta_norm*theta_norm) - (1+std::cos(theta_norm)) / (2*theta_norm*std::sin(theta_norm)) ) * skew * skew;
}

static Mat3r Exp_so3(const Vec3r& vec)
{
    const Mat3r skew = Skew3(vec);
    Real mag = vec.norm();

    if (mag < Real(1e-8))
        return Mat3r::Identity() + skew;
    
    return Mat3r::Identity() + std::sin(mag) / mag * skew + (1 - std::cos(mag)) / (mag * mag) * skew * skew;
}

static Vec3r Log_SO3(const Mat3r& mat)
{
    // std::cout << "\n===Log_SO3===" << std::endl;
    // std::cout << "  mat:\n" << mat << std::endl;
    // std::cout << "  mat.trace()-3: " << mat.trace()-3 << std::endl;
    Real theta = std::acos( std::min(0.5 * mat.trace() - 0.5, Real(1.0)));  // make sure 1/2 tr(mat) - 1/2 is not >1, will get NaNs. This may happen due to numerical drift
    // std::cout << "  theta: " << theta << std::endl;

    if (std::abs(theta) < Real(1e-14))
    {
        return Vec3r::Zero();
    }

    const Vec3r skew_vec3 = Vee3(mat - mat.transpose());
    // std::cout << "  skew_vec3: " << skew_vec3 << std::endl;

    if (std::abs(mat.trace()) < Real(1e-8))
    {
        return 0.5 * (1 + theta*theta/6.0 + 7*theta*theta*theta*theta/360.0) * skew_vec3;
    }

    // std::cout << " 2*std::sin(theta): " << 2*std::sin(theta) << std::endl;
    return theta / ( 2*std::sin(theta)) * skew_vec3;
}

static Vec3r Minus_SO3(const Mat3r& mat1, const Mat3r& mat2)
{
    return Log_SO3(mat2.transpose() * mat1);
}

static Mat3r Plus_SO3(const Mat3r& SO3_mat, const Vec3r& so3_vec)
{
    return SO3_mat * Exp_so3(so3_vec);
}


/** === Quaternions === */
static Quaternion Exp_s3(const Vec3r& vec)
{
    Real theta = vec.norm();
    Quaternion q;

    // small angle approximation
    if (theta < 1e-6)
    {
        q.w() = 1.0;
        q.x() = 0.5 * vec.x();
        q.y() = 0.5 * vec.y();
        q.z() = 0.5 * vec.z();
        q.normalize();
        return q;
    }

    Real half_theta = 0.5*theta;
    Real sin_half_theta = std::sin(half_theta);
    Real scale = sin_half_theta / theta;
    q.w() = std::cos(half_theta);
    q.x() = scale * vec.x();
    q.y() = scale * vec.y();
    q.z() = scale * vec.z();

    return q;
}

static Vec3r Log_S3(const Quaternion& q)
{
    Vec3r v(q.x(), q.y(), q.z());
    Real v_mag = v.norm();
    if (v_mag < 1e-6)
    {
        return 2*v;
    }

    Real theta = 2*std::atan2(v_mag, q.w());
    return theta/v_mag * v;
    
}

static Vec3r Minus_S3(const Quaternion& q1, const Quaternion& q2)
{
    return Log_S3(q2.conjugate() * q1);
}

static Quaternion Plus_S3(const Quaternion& q, const Vec3r& rot_vec)
{
    return q * Exp_s3(rot_vec);
}

static Mat3r RotMatFromXYZEulerAngles(const Vec3r& euler_xyz)
{
    const Real x = euler_xyz(0) * M_PI / 180.0;
    const Real y = euler_xyz(1) * M_PI / 180.0;
    const Real z = euler_xyz(2) * M_PI / 180.0;

    // using the "123" convention: rotate first about x axis, then about y, then about z
    Mat3r rot_mat;
    rot_mat(0,0) = std::cos(y) * std::cos(z);
    rot_mat(0,1) = std::sin(x)*std::sin(y)*std::cos(z) - std::cos(x)*std::sin(z);
    rot_mat(0,2) = std::cos(x)*std::sin(y)*std::cos(z) + std::sin(x)*std::sin(z);

    rot_mat(1,0) = std::cos(y)*std::sin(z);
    rot_mat(1,1) = std::sin(x)*std::sin(y)*std::sin(z) + std::cos(x)*std::cos(z);
    rot_mat(1,2) = std::cos(x)*std::sin(y)*std::sin(z) - std::sin(x)*std::cos(z);

    rot_mat(2,0) = -std::sin(y);
    rot_mat(2,1) = std::sin(x)*std::cos(y);
    rot_mat(2,2) = std::cos(x)*std::cos(y);

    return rot_mat;
}

static Vec3r XYZEulerAnglesFromRotMat(const Mat3r& R)
{
    Real theta_y = std::asin(-R(2,0));
    Real theta_x, theta_z;
    if (std::abs(R(2,0)) < 1)
    {
        theta_x = std::atan2(R(2,1), R(2,2));
        theta_z = std::atan2(R(1,0), R(0,0));
    }
    else
    {
        theta_x = 0;
        if (R(2,0) == -1)
        {
            theta_y = M_PI/2;
            theta_z = std::atan2(-R(0,1), R(0,2));
        }
        else
        {
            theta_y = -M_PI/2;
            theta_z = std::atan2(R(0,1), -R(0,2));
        }
    }

    return 180/M_PI * Vec3r(theta_x, theta_y, theta_z);
}

static Quaternion QuaternionFromXYZEulerAngles(const Vec3r& eul_xyz)
{
    Real x = eul_xyz[0] * M_PI / 180.0;
    Real y = eul_xyz[1] * M_PI / 180.0;
    Real z = eul_xyz[2] * M_PI / 180.0;

    Real cx = std::cos(x * 0.5);
    Real sx = std::sin(x * 0.5);

    Real cy = std::cos(y * 0.5);
    Real sy = std::sin(y * 0.5);

    Real cz = std::cos(z * 0.5);
    Real sz = std::sin(z * 0.5);

    Quaternion qx(sx, 0, 0, cx);
    Quaternion qy(0, sy, 0, cy);
    Quaternion qz(0, 0, sz, cz);

    return qz*qy*qx;
}

/** Given a unit normal vector, completes an orthonormal basis.
 * Uses Eigen's unitOrthogonal().
 */
static void completeOrthonormalBasisGivenNormal(const Vec3r& n, Vec3r& t, Vec3r& b)
{
    t = n.unitOrthogonal();
    b = n.cross(t);
}

/** Projects a point p onto the line segment defined by ab.
 * Returns the interpolation factor - i.e. if in [0,1] the projected point is between a and b.
 * To get the projected point: p_proj = a + projectPointOnLine(p, a, b) * (b-a);
 */
static Real projectPointOntoLine(const Vec3r& p, const Vec3r& a, const Vec3r& b)
{
    return (p-a).dot(b-a) / (b-a).squaredNorm();
}

/** Find closest points on two line segments defined by (p1, p2) and (p3, p4) */
static std::pair<Real, Real> findClosestPointsOnLineSegments(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4)
{
    const Vec3r d1 = p2 - p1;
    const Vec3r d2 = p4 - p3;
    const Vec3r w = p1 - p3;

    Real a = d1.dot(d1);
    Real b = d1.dot(d2);
    Real c = d2.dot(d2);
    Real d = w.dot(d1);
    Real e = w.dot(d2);
    
    Real den = a*c - b*b;

    Real beta1, beta2;
    if (den < 1e-8) // line segments are roughly parallel
    {
        Real best_dist_sq = std::numeric_limits<Real>::max();

        // project p1 onto (p3,p4)
        Real beta2_1 = std::clamp(projectPointOntoLine(p1, p3, p4), 0.0, 1.0);
        Real dist_sq1 = ((1-beta2_1)*p3 + beta2_1*p4 - p1).squaredNorm();
        best_dist_sq = dist_sq1;
        beta1 = 0.0;
        beta2 = beta2_1;

        // project p2 onto (p3, p4)
        Real beta2_2 = std::clamp(projectPointOntoLine(p2, p3, p4), 0.0, 1.0);
        Real dist_sq2 = ((1-beta2_2)*p3 + beta2_2*p4 - p2).squaredNorm();
        if (dist_sq2 < best_dist_sq)
        {
            best_dist_sq = dist_sq2;
            beta1 = 1.0;
            beta2 = beta2_2;
        }

        // project p3 onto (p1, p2)
        Real beta1_1 = std::clamp(projectPointOntoLine(p3, p1, p2), 0.0, 1.0);
        Real dist_sq3 = ((1-beta1_1)*p1 + beta1_1*p2 - p3).squaredNorm();
        if (dist_sq3 < best_dist_sq)
        {
            best_dist_sq = dist_sq3;
            beta1 = beta1_1;
            beta2 = 0.0;
        }

        // project p4 onto (p1, p2)
        Real beta1_2 = std::clamp(projectPointOntoLine(p4, p1, p2), 0.0, 1.0);
        Real dist_sq4 = ((1-beta1_2)*p1 + beta1_2*p2 - p4).squaredNorm();
        if (dist_sq4 < best_dist_sq)
        {
            best_dist_sq = dist_sq4;
            beta1 = beta1_2;
            beta2 = 1.0;
        }
    }
    else
    {
        beta1 = (b*e - c*d) / den;
        beta2 = (a*e - b*d) / den;

        if (beta1 < 0 || beta1 > 1)
        {
            // if beta1 was clamped, re-project to find beta2
            beta1 = std::clamp(beta1, 0.0, 1.0);
            beta2 = std::clamp(projectPointOntoLine((1-beta1)*p1 + beta1*p2, p3, p4), 0.0, 1.0);
            
        }
        else if (beta2 < 0 || beta2 > 1)
        {
            // if beta2 was clamped, re-project to find beta1
            beta2 = std::clamp(beta2, 0.0, 1.0);
            beta1 = std::clamp(projectPointOntoLine((1-beta2)*p3 + beta2*p4, p1, p2), 0.0, 1.0);
        }
    }

    return std::make_pair(beta1, beta2);
}


/** Geometry subroutines */

inline static Vec3r barycentricCoordinates(const Vec3r& p, const Vec3r& a, const Vec3r& b, const Vec3r& c)
{
    // from https://ceng2.ktu.edu.tr/~cakir/files/grafikler/Texture_Mapping.pdf
    const Vec3r v0 = b - a;
    const Vec3r v1 = c - a;
    const Vec3r v2 = p - a;
    const Real d00 = v0.dot(v0);
    const Real d01 = v0.dot(v1);
    const Real d11 = v1.dot(v1);
    const Real d20 = v2.dot(v0);
    const Real d21 = v2.dot(v1);
    const Real denom = d00*d11 - d01*d01;
    const Real v = (d11*d20 - d01*d21) / denom;
    const Real w = (d00*d21 - d01*d20) / denom;
    const Real u = 1 - v - w;

    return Vec3r(u, v, w);
}

/** Closest point on triangle to query point.
 * @param p the query point
 * @param a,b,c triangle vertices 
 * 
 * adapted from: https://github.com/RenderKit/embree/blob/master/tutorials/common/math/closest_point.h
 */
inline static Vec3r closestPoint_PointTriangle(const Vec3r& p, const Vec3r& a, const Vec3r& b, const Vec3r& c)
{
    const Vec3r ab = b - a;
    const Vec3r ac = c - a;
    const Vec3r ap = p - a;

    const Real d1 = ab.dot(ap);
    const Real d2 = ac.dot(ap);
    if (d1 <= 0.f && d2 <= 0.f)
    {
        return a;
    }

    const Vec3r bp = p - b;
    const Real d3 = ab.dot(bp);
    const Real d4 = ac.dot(bp);
    if (d3 >= 0.f && d4 <= d3)
    {
        return b;
    }

    const Vec3r cp = p - c;
    const Real d5 = ab.dot(cp);
    const Real d6 = ac.dot(cp);
    if (d6 >= 0.f && d5 <= d6)
    {
        return c;
    }

    const Real vc = d1 * d4 - d3 * d2;
    if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
    {
        const Real v = d1 / (d1 - d3);
        return a + v * ab;
    }
    
    const Real vb = d5 * d2 - d1 * d6;
    if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
    {
        const Real v = d2 / (d2 - d6);
        return a + v * ac;
    }
    
    const Real va = d3 * d6 - d5 * d4;
    if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
    {
        const Real v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + v * (c - b);
    }

    const Real denom = 1.f / (va + vb + vc);
    const Real v = vb * denom;
    const Real w = vc * denom;
    return a + v * ab + w * ac;
}

/** Closest points between two line segments
 * Adapted from Christer Ericson's ClosestPtSegmentSegment() https://ceng2.ktu.edu.tr/~cakir/files/grafikler/rtcd.pdf
 * @param p1,q1 : line segment 1
 * @param p2,q2 : line segment 2
 * @param s,t (output) the interpolation parameters for the closest points on segment 1 and segment 2, respectively
 */
inline static void closestPoint_SegmentSegment(const Vec3r& p1, const Vec3r& q1, const Vec3r& p2, const Vec3r& q2, Real& s, Real& t)
{
    Vec3r d1 = q1 - p1; // Direction vector of segment S1
    Vec3r d2 = q2 - p2; // Direction vector of segment S2
    Vec3r r = p1 - p2;
    Real a = d1.squaredNorm(); // Squared length of segment S1, always nonnegative
    Real e = d2.squaredNorm(); // Squared length of segment S2, always nonnegative
    Real f = d2.dot(r);

    // Check if either or both segments degenerate into points
    Real EPSILON = 1e-8;
    if (a <= EPSILON && e <= EPSILON) {
        // Both segments degenerate into points
        s = t = 0.0;
        return;
    }
    if (a <= EPSILON) {
        // First segment degenerates into a point
        s = 0.0;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = std::clamp(t, Real(0.0), Real(1.0));
    } else {
        Real c = d1.dot(r);
        if (e <= EPSILON) {
            // Second segment degenerates into a point
            t = 0.0f;
            s = std::clamp(-c / a, Real(0.0), Real(1.0)); // t = 0 => s = (b*t - c) / a = -c / a
        } else {
            // The general nondegenerate case starts here
            Real b = d1.dot(d2);
            Real denom = a*e-b*b; // Always nonnegative
            // If segments not parallel, compute closest point on L1 to L2 and
            // clamp to segment S1. Else pick arbitrary s (here 0)
            if (denom != 0.0) {
                s = std::clamp((b*f - c*e) / denom, Real(0.0), Real(1.0));
            } else s = 0.0;

            // Compute point on L2 closest to S1(s) using
            // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
            t = (b*s + f) / e;
            // If t in [0,1] done. Else clamp t, recompute s for the new value
            // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a
            // and clamp s to [0, 1]
            if (t < Real(0.0)) {
                t = 0.0;
                s = std::clamp(-c / a, Real(0.0), Real(1.0));
            } else if (t > 1.0) {
                t = 1.0;
                s = std::clamp((b - c) / a, Real(0.0), Real(1.0));
            }
        }
    }
}
};