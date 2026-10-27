#pragma once

#include "energy/collision/CollisionConstraintEnergyPool.hpp"
#include "energy/collision/TriangleRodCollisionEnergyDefs.hpp"

namespace Energy
{

struct TriangleRodCollisionEnergyPool : CollisionConstraintEnergyPool<TriangleRodCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 5;
    static constexpr EnergyType Type = EnergyType::TRIANGLE_ROD_COLLISION;
    static constexpr DynamicEnergyType DynamicType = DynamicEnergyType::TRIANGLE_ROD_COLLISION;
    using SolverType = TriangleRodCollisionEnergySolver;
    using ViewType = TriangleRodCollisionEnergyPoolView;

    explicit TriangleRodCollisionEnergyPool(unsigned capacity)
        : CollisionConstraintEnergyPool(capacity)
    {
    }

    TriangleRodCollisionEnergyPoolView view()
    {
        TriangleRodCollisionEnergyPoolView v;
        CollisionConstraintEnergyPool::view(v);
        return v;
    }

    /** Add an energy
     * @param p_idx1, p_idx2, p_idx3 the particle indices of the triangle face
     * @param s_idx1, s_idx2 the particle indices of the rod segment
     * @param sdf_params a pointer to the SDF parameters for the rigid body
     */
    unsigned addEnergy(
        unsigned p_idx1, 
        unsigned p_idx2, 
        unsigned p_idx3, 
        unsigned s_idx1,
        unsigned s_idx2, 
        const Vec3r& normal,
        const Vec3r& tangent,
        const Vec3r& binormal,
        const Vec3r& barys,
        Real s,
        const Vec3r& cp_rod_local,
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
        data[slot].particle_indices[3] = s_idx1;
        data[slot].particle_indices[4] = s_idx2;
        data[slot].normal = normal;
        data[slot].tangent = tangent;
        data[slot].binormal = binormal;
        data[slot].barys = barys;
        data[slot].s = s;
        data[slot].cp_rod_local = cp_rod_local;

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