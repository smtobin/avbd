#pragma once

#include "energy/CollisionConstraintEnergyPool.hpp"
#include "collision/SDF.hpp"

namespace Energy
{

struct TriangleRodCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec5u particle_indices;     // particle indices - first 3 indices are the triangle vertices, last 2 indices make up the rod segment
    Vec3r barys;    // barycentric coordinates of the contact point on the face
    Real s;  // location of contact point along rod segment (in the interval [0,1])
    Vec3r cp_rod_local; // local contact point in the rod interpolated frame
};

struct TriangleRodCollisionEnergyPool : CollisionConstraintEnergyPool<TriangleRodCollisionEnergyInfo>
{
    static constexpr int NumParticlesPerEnergy = 5;
    static constexpr EnergyType Type = EnergyType::TRIANGLE_ROD_COLLISION;
    static constexpr DynamicEnergyType DynamicType = DynamicEnergyType::TRIANGLE_ROD_COLLISION;
    using SolverType = TriangleRodCollisionEnergySolver;

    explicit TriangleRodCollisionEnergyPool(unsigned capacity)
        : CollisionConstraintEnergyPool(capacity, 1e2)
    {
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
        Real mu_s,
        Real mu_k
    )
    {
        // parent will call allocSlot()
        unsigned slot = CollisionConstraintEnergyPool::addEnergy(mu_s, mu_k);

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