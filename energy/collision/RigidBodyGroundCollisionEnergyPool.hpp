#pragma once

#include "energy/collision/CollisionConstraintEnergyPool.hpp"
#include "energy/collision/RigidBodyGroundCollisionEnergyDefs.hpp"

namespace Energy
{

struct RigidBodyGroundCollisionEnergyPool : CollisionConstraintEnergyPool<RigidBodyGroundCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 1; // number of particles per constraint
    static constexpr EnergyType Type = EnergyType::RIGID_BODY_GROUND_COLLISION; // type of energy in the EnergyType enum
    static constexpr StaticEnergyType StaticType = StaticEnergyType::RIGID_BODY_GROUND_COLLISION; // type of energy in the StaticEnergyType enum
    using SolverType = RigidBodyGroundCollisionEnergySolver;     // solver class type

    explicit RigidBodyGroundCollisionEnergyPool(unsigned capacity)
        : CollisionConstraintEnergyPool(capacity)
    {

    }

    /** Returns a non-owning view into the memory pool. */
    RigidBodyGroundCollisionEnergyPoolView view() const
    {
        RigidBodyGroundCollisionEnergyPoolView v;
        CollisionConstraintEnergyPool::view(v);
        return v;
    }

    /** Add an energy
     * @param particle_index : the index of the particle in the particle pool
     * @param sdf_index : the index of the rigid body in the sim collision pool's SDF pool
     * @param cp_rb_local : the initial body-frame vector to the contact point on the rigid body
     * @returns the index of the new energy in the pool
     */
    unsigned addEnergy(
        unsigned particle_index, 
        const Collision::CollisionShapeParams& shape_params,
        const Vec3r& cp_rb_local, 
        Real k_start, 
        Real mu_s, 
        Real mu_k
    )
    {
        // parent will call allocSlot()
        unsigned slot = CollisionConstraintEnergyPool::addEnergy(k_start, mu_s, mu_k);

        // particle_indices[slot][0] = particle_index;
        data[slot].particle_indices[0] = particle_index;
        data[slot].shape_params = shape_params;
        data[slot].cp_rb_local = cp_rb_local;
        data[slot].cp_ground = Vec3r::Zero();
        data[slot].normal = Vec3r(0,1,0);
        data[slot].tangent = Vec3r(1,0,0);
        data[slot].binormal = Vec3r(0,0,1);

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