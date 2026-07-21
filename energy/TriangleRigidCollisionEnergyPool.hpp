#pragma once

#include "energy/HardConstraintEnergyPool.hpp"
#include "collision/SDF.hpp"

namespace Energy
{

struct TriangleRigidCollisionEnergyInfo : HardConstraintEnergyInfo
{
    Vec3u particle_indices;
    Vec1u oriented_particle_indices;
    Collision::SDFShapeParams* sdf_params;
};

struct TriangleRigidCollisionEnergyPool : HardConstraintEnergyPool<TriangleRigidCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 3;
    static constexpr EnergyType Type = EnergyType::TRIANGLE_RIGID_COLLISION;
    using SolverType = TriangleRigidCollisionEnergySolver;

    explicit TriangleRigidCollisionEnergyPool(unsigned capacity)
        : HardConstraintEnergyPool(capacity, 1e2, 0)
    {
    }

    /** Add an energy
     * @param p_idx1, p_idx2, p_idx3 the particle indices of the triangle face
     * @param op_idx the particle index of the rigid body
     * @param sdf_params a pointer to the SDF parameters for the rigid body
     */
    unsigned addEnergy(unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, unsigned op_idx, Collision::SDFShapeParams* sdf_params)
    {
        // parent will call allocSlot()
        unsigned slot = HardConstraintEnergyPool::addEnergy();

        data[slot].particle_indices[0] = p_idx1;
        data[slot].particle_indices[1] = p_idx2;
        data[slot].particle_indices[2] = p_idx3;
        data[slot].oriented_particle_indices[0] = op_idx;
        data[slot].sdf_params = sdf_params;

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