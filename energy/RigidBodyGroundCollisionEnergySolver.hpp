#pragma once

#include "energy/RigidBodyGroundCollisionEnergyPool.hpp"
#include "energy/CollisionConstraintEnergySolver.hpp"

namespace Energy
{

struct RigidBodyGroundCollisionConstraintSolver
{
    /** Public typedefs */
    using PoolType = RigidBodyGroundCollisionEnergyPool;
    
    static void evaluateConstraint(
        unsigned c_idx,
        const PoolType& energies,
        ParticlePool& particles,
        Real& C_n, Real& C_t, Real& C_b
    )
    {
        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3R& cp_ground = energies.data[c_idx].cp_ground;
        Vec3r cp_rb = particles.positions[p_idx] + particles.rotation(p_idx) * cp_rb_local;
        Vec3r diff = cp_rb - cp_ground;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);
    }

    static void constraintGradientHessian(
        unsigned c_idx,
        const PoolType& energies,
        ParticlePool& particles,
        unsigned /* local_idx */,
        Real& C_n, Real& C_t, Real& C_b,
        Vec6r<DOF>& C_grad_n, Vec3r_or_Vec6r<DOF>& C_grad_t, Vec3r_or_Vec6r<DOF>& C_grad_b,
        Mat6r<DOF>& C_hess_n, Mat3r_or_Mat6r<DOF>& C_hess_t, Mat3r_or_Mat6r<DOF>& C_hess_b
    )
    {
        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3r& cp_ground = energies.data[c_idx].cp_ground;
        const Quaternion& rb_rot = particles.rotation(p_idx);
        Vec3r cp_rb = particles.positions[p_idx] + particles.rotation(p_idx) * cp_rb_local;
        Vec3r diff = cp_rb - cp_ground;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);
        
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
};

// using GroundCollisionEnergySolver = 
//     HardConstraintEnergySolver<GroundCollisionEnergyPool, GroundCollisionConstraintSolver>;

struct RigidBodyGroundCollisionEnergySolver
    : CollisionConstraintEnergySolver<RigidBodyGroundCollisionEnergyPool, RigidBodyGroundCollisionConstraintSolver>
{
    using CollisionConstraintEnergySolver::CollisionConstraintEnergySolver;
    static constexpr bool SupportsPositional = false;
    static constexpr bool SupportsOriented = true;

    /** Updates the stiffness and Lagrange multiplier after a full time step.
     * Implements equation (19) from the AVBD paper.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     */
    static void updateAfterTimeStep(
        unsigned c_idx,
        RigidBodyGroundCollisionEnergyPool& energies,
        ParticlePool& particles
    )
    {
        // update the contact point
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        energies.data[c_idx].cp_ground = particles.positions[p_idx];

        // update the local contact point, when necessary
        // (for spheres and capsules)
        Collision::SDFShapeParams* sdf_params = energies.data[c_idx].sdf_params;
        if (sdf_params->type == Collision::SDFType::Sphere)
        {
            energies.data[c_idx].cp_rb_local = particles.rotation(p_idx).conjugate() * Vec3r(0, -sdf_params->sphere.radius, 0);
        }

        CollisionConstraintEnergySolver::updateAfterTimeStep(c_idx, energies, particles);
    }
};

} // namespace Energy