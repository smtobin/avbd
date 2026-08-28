#pragma once

#include "energy/collision/RigidBodyGroundCollisionEnergyPool.hpp"
#include "energy/collision/CollisionConstraintEnergySolver.hpp"

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
        const Collision::SDFShapeParams& shape_params = energies.data[c_idx].shape_params;
        unsigned p_idx = energies.data[c_idx].particle_indices[0];
        
        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3r& cp_ground = energies.data[c_idx].cp_ground;
        Vec3r cp_rb = particles.positions[p_idx] + particles.rotation(p_idx) * cp_rb_local;
        Vec3r diff = cp_ground - cp_rb;
        // C_n = n.dot(diff);
        C_t = t.dot(diff);
        C_b = b.dot(diff);

        // if shape is round (e.g. sphere, rod), always use the exact contact point
        // (instead of the one computed at the beginning of the time step)
        // this will always be a distance of radius in the -y direction
        // this accommodates high rolling and prevents weird behavior with rod-ground collisions
        if (shape_params.type == Collision::SDFType::Sphere)
        {
            Vec3r cp_global = particles.positions[p_idx] - Vec3r(0, shape_params.sphere.radius, 0);
            C_n = n.dot(cp_ground - cp_global);
        }
        else if (shape_params.type == Collision::SDFType::Rod)
        {
            Vec3r cp_global = particles.positions[p_idx] - Vec3r(0, shape_params.rod.radius, 0);
            C_n = n.dot(cp_ground - cp_global);
        }
        else
        {
            // default, just use the contact point
            C_n = n.dot(diff);
        }
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
        const Collision::SDFShapeParams& shape_params = energies.data[c_idx].shape_params;
        unsigned p_idx = energies.data[c_idx].particle_indices[0];

        const Vec3r& cp_rb_local = energies.data[c_idx].cp_rb_local;
        const Vec3r& cp_ground = energies.data[c_idx].cp_ground;
        const Quaternion& rb_rot = particles.rotation(p_idx);
        Vec3r cp_rb = particles.positions[p_idx] + particles.rotation(p_idx) * cp_rb_local;

        /** TEST! */
        // cp_rb = particles.positions[p_idx] + Vec3r(0, -0.5, 0);


        Vec3r diff = cp_ground - cp_rb;
        C_t = t.dot(diff);
        C_b = b.dot(diff);

        // precomputations for gradients
        Mat3r R = rb_rot.toRotationMatrix();
        Mat3r skew_cp = Math::Skew3(cp_rb_local);
        Mat3r R_rloc = R * skew_cp;

        // if shape is round (e.g. sphere, rod), always use the exact contact point
        // (instead of the one computed at the beginning of the time step)
        // this will always be a distance of radius in the -y direction
        // this accommodates high rolling and prevents weird behavior with rod-ground collisions
        if (shape_params.type == Collision::SDFType::Sphere)
        {
            Vec3r cp_global = particles.positions[p_idx] - Vec3r(0, shape_params.sphere.radius, 0);
            C_n = n.dot(cp_ground - cp_global);
            C_grad_n.template block<3,1>(0,0) = -n;
            C_grad_n.template block<3,1>(3,0) = Vec3r::Zero();
            C_hess_n = Mat6r::Zero();
        }
        else if (shape_params.type == Collision::SDFType::Rod)
        {
            Vec3r cp_global = particles.positions[p_idx] - Vec3r(0, shape_params.rod.radius, 0);
            C_n = n.dot(cp_ground - cp_global);
            C_grad_n.template block<3,1>(0,0) = -n;
            C_grad_n.template block<3,1>(3,0) = Vec3r::Zero();
            C_hess_n = Mat6r::Zero();
        }
        else
        {
            // default, just use the contact point
            C_n = n.dot(diff);
            C_grad_n.template block<3,1>(0,0) = -n;
            C_grad_n.template block<3,1>(3,0) = n.transpose() * R_rloc;
            C_hess_n = Mat6r::Zero();
            C_hess_n.template block<3,3>(3,3) = -R.transpose() * Math::Skew3(n) * R * skew_cp;
            C_hess_n.template block<3,3>(3,3) = 0.5*(C_hess_n.template block<3,3>(3,3) + C_hess_n.template block<3,3>(3,3).transpose());
        }

        // std::cout << "C_n: " << C_n << "  C_t: " << C_t << "  C_b: " << C_b << std::endl;
        // std::cout << "cp_rb_local: " << cp_rb_local.transpose() << "\nreal cp_rb_local: " << (rb_rot.inverse() * Vec3r(0, -0.5, 0)).transpose() << std::endl;
        
        
        // gradients
        C_grad_t.template block<3,1>(0,0) = -t;
        C_grad_t.template block<3,1>(3,0) = t.transpose() * R_rloc;
        C_grad_b.template block<3,1>(0,0) = -b;
        C_grad_b.template block<3,1>(3,0) = b.transpose() * R_rloc;

        // Hessians
        // grad = (Skew(cp_rb_local) * R^T * n)^T
        // ==> (Skew(cp_rb_local) * -R^T * skew(n) * -R)^T = R^T * skew(n) * R * skew(cp_rb_local)
        C_hess_t = C_hess_b = Mat6r::Zero();
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
        const Collision::SDFShapeParams& sdf_params = energies.data[c_idx].shape_params;
        if (sdf_params.type == Collision::SDFType::Sphere)
        {
            energies.data[c_idx].cp_rb_local = particles.rotation(p_idx).conjugate() * Vec3r(0, -sdf_params.sphere.radius, 0);
        }
        else if (sdf_params.type == Collision::SDFType::Rod)
        {
            energies.data[c_idx].cp_rb_local = particles.rotation(p_idx).conjugate() * Vec3r(0, -sdf_params.rod.radius, 0);
        }
    }
};

} // namespace Energy