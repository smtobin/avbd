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
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec3r& t = energies.data[c_idx].tangent;
        const Vec3r& b = energies.data[c_idx].binormal;
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        
        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3r& cp_ground = energies.data[c_idx].cp_ground;
        Vec3r cp_rb = particles.positions[p_idx] + particles.rotation(p_idx) * cp_rb_local;
        Vec3r diff = cp_ground - cp_rb;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);
    }

    template <int>
    static void constraintGradientHessian(
        unsigned c_idx,
        const PoolType& energies,
        ParticlePool& particles,
        unsigned /* local_idx */,
        Real& C_n, Real& C_t, Real& C_b,
        Vec6r& C_grad_n, Vec6r& C_grad_t, Vec6r& C_grad_b,
        Mat6r& C_hess_n, Mat6r& C_hess_t, Mat6r& C_hess_b
    )
    {
        const Vec3r& n = energies.data[c_idx].normal;
        const Vec3r& t = energies.data[c_idx].tangent;
        const Vec3r& b = energies.data[c_idx].binormal;
        unsigned p_idx = energies.data[c_idx].particle_indices[0];

        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3r& cp_ground = energies.data[c_idx].cp_ground;
        const Quaternion& rb_rot = particles.rotation(p_idx);
        Vec3r cp_rb = particles.positions[p_idx] + particles.rotation(p_idx) * cp_rb_local;
        Vec3r diff = cp_ground - cp_rb;
        C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);
        
        Mat3r R = rb_rot.toRotationMatrix();
        Mat3r skew_cp = Math::Skew3(cp_rb_local);
        Mat3r R_rloc = R * skew_cp;
        // gradients
        C_grad_n.template block<3,1>(0,0) = -n;
        C_grad_n.template block<3,1>(3,0) = n.transpose() * R_rloc;
        C_grad_t.template block<3,1>(0,0) = -t;
        C_grad_t.template block<3,1>(3,0) = t.transpose() * R_rloc;
        C_grad_b.template block<3,1>(0,0) = -b;
        C_grad_b.template block<3,1>(3,0) = b.transpose() * R_rloc;
        C_grad_t = C_grad_b = Vec6r::Zero();

        // Hessians
        // grad = (Skew(cp_rb_local) * R^T * n)^T
        // ==> (Skew(cp_rb_local) * -R^T * skew(n) * -R)^T = R^T * skew(n) * R * skew(cp_rb_local)
        C_hess_n = C_hess_t = C_hess_b = Mat6r::Zero();
        C_hess_n.template block<3,3>(3,3) = -R.transpose() * Math::Skew3(n) * R * skew_cp;
        C_hess_n.template block<3,3>(3,3) = 0.5*(C_hess_n.template block<3,3>(3,3) + C_hess_n.template block<3,3>(3,3).transpose());
        C_hess_t.template block<3,3>(3,3) = -R.transpose() * Math::Skew3(t) * R * skew_cp;
        C_hess_t.template block<3,3>(3,3) = 0.5*(C_hess_t.template block<3,3>(3,3) + C_hess_t.template block<3,3>(3,3).transpose());
        C_hess_b.template block<3,3>(3,3) = -R.transpose() * Math::Skew3(b) * R * skew_cp;
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

    /** Updates the contact point for this energy.
     * Since this energy is static (created upon simulation initialization), it is not updated by collision detection.
     * So, we must update the contact point here - this is called at the start of each time step.
     * @param c_idx : the constraint index
     * @param energies : the memory pool for the energy
     * @param particles : the particle pool
     * @param sdf_pool : the SDF primitive pool
     */
    static void updateContactPoints(
        unsigned c_idx,
        RigidBodyGroundCollisionEnergyPool& energies,
        ParticlePool& particles,
        Collision::SDFPrimitivePool& sdf_pool
    )
    {
        // update the contact point
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        energies.data[c_idx].cp_ground[0] = particles.positions[p_idx][0];
        energies.data[c_idx].cp_ground[2] = particles.positions[p_idx][2];

        // update the local contact point, when necessary
        // (for spheres and capsules)
        unsigned sdf_idx = energies.data[c_idx].sdf_index;
        const Collision::SDFShapeParams& sdf_params = sdf_pool.params[sdf_idx];
        if (sdf_params.type == Collision::SDFType::Sphere)
        {
            energies.data[c_idx].cp_rb_local = particles.rotation(p_idx).conjugate() * Vec3r(0, -sdf_params.sphere.radius, 0);
        }
    }
};

} // namespace Energy