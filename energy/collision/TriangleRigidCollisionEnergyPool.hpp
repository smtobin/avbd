#pragma once

// #include "energy/HardConstraintEnergyPool.hpp"
#include "energy/collision/CollisionConstraintEnergyPool.hpp"
#include "collision/SDF.hpp"

namespace Energy
{

struct TriangleRigidCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec4u particle_indices;     // particle indices - first 3 indices are the triangle vertices, last index is the rigid body
    Collision::CollisionShapeParams* sdf_params;  // pointer to SDF parameters - will be used to reevaluate the constraint
    Vec3r barys;    // barycentric coordinates of the contact point on the face
    Vec3r cp_rb_local;  // contact point on the rigid body, in the rigid body's local frame
};

struct TriangleRigidCollisionEnergyPool : CollisionConstraintEnergyPool<TriangleRigidCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 4;
    static constexpr EnergyType Type = EnergyType::TRIANGLE_RIGID_COLLISION;
    static constexpr DynamicEnergyType DynamicType = DynamicEnergyType::TRIANGLE_RIGID_COLLISION;
    using SolverType = TriangleRigidCollisionEnergySolver;

    explicit TriangleRigidCollisionEnergyPool(unsigned capacity)
        : CollisionConstraintEnergyPool(capacity)
    {
    }

    /** Add an energy
     * @param p_idx1, p_idx2, p_idx3 the particle indices of the triangle face
     * @param op_idx the particle index of the rigid body
     * @param sdf_params a pointer to the SDF parameters for the rigid body
     */
    unsigned addEnergy(
        unsigned p_idx1, 
        unsigned p_idx2, 
        unsigned p_idx3, 
        unsigned op_idx, 
        Collision::CollisionShapeParams* sdf_params,
        const Vec3r& normal,
        const Vec3r& tangent,
        const Vec3r& binormal,
        const Vec3r& barys,
        const Vec3r& cp_rb_local,
        Real k_start, 
        Real mu_s,
        Real mu_k
    )
    {
        // parent will call allocSlot()
        unsigned slot = CollisionConstraintEnergyPool::addEnergy(k_start, mu_s, mu_k);

        data[slot].particle_indices[0] = p_idx1;
        data[slot].particle_indices[1] = p_idx2;
        data[slot].particle_indices[2] = p_idx3;
        data[slot].particle_indices[3] = op_idx;
        data[slot].sdf_params = sdf_params;
        data[slot].normal = normal;
        data[slot].tangent = tangent;
        data[slot].binormal = binormal;
        data[slot].barys = barys;
        data[slot].cp_rb_local = cp_rb_local;

        return slot;
    }

    /** Remove an energy
     * @param slot : the index of the energy in the pool to remove
     */
    void removeEnergy(unsigned slot)
    {
        freeSlot(slot);
    }
};

} // namespace Energy