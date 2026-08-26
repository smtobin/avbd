#pragma once

#include "energy/NeoHookeanEnergyPool.hpp"
#include "energy/CosseratRodEnergyPool.hpp"
#include "energy/GroundCollisionEnergyPool.hpp"
#include "energy/RigidBodyGroundCollisionEnergyPool.hpp"
#include "energy/TriangleRigidCollisionEnergyPool.hpp"
#include "energy/TriangleRodCollisionEnergyPool.hpp"
#include "energy/OneSidedFixedJointEnergyPool.hpp"

namespace Energy
{

/** Storage of all the energies in the sim.
 * Stores the individual energy memory pools.
 */
struct EnergyRegistry
{
    // elastic
    NeoHookeanEnergyPool neo_hookean;
    CosseratRodEnergyPool cosserat_rod;

    // collision
    GroundCollisionEnergyPool ground_collision;
    RigidBodyGroundCollisionEnergyPool rigid_body_ground_collision;
    TriangleRigidCollisionEnergyPool triangle_rigid_collision;
    TriangleRodCollisionEnergyPool triangle_rod_collision;

    // joints
    OneSidedFixedJointEnergyPool one_sided_fixed_joint;

    EnergyRegistry(unsigned capacity)
        : neo_hookean(capacity)
        , cosserat_rod(capacity)
        , ground_collision(capacity)
        , rigid_body_ground_collision(capacity)
        , triangle_rigid_collision(capacity)
        , triangle_rod_collision(capacity)
        , one_sided_fixed_joint(capacity)
    {}

    /** Statically-typed getter for a specific energy type */
    template<EnergyType E>
    auto& get()
    {
        if constexpr (E == EnergyType::NEO_HOOKEAN)
            return neo_hookean;
        else if constexpr (E == EnergyType::COSSERAT_ROD)
            return cosserat_rod;
        else if constexpr (E == EnergyType::GROUND_COLLISION)
            return ground_collision;
        else if constexpr (E == EnergyType::RIGID_BODY_GROUND_COLLISION)
            return rigid_body_ground_collision;
        else if constexpr (E == EnergyType::TRIANGLE_RIGID_COLLISION)
            return triangle_rigid_collision;
        else if constexpr (E == EnergyType::TRIANGLE_ROD_COLLISION)
            return triangle_rod_collision;
        else if constexpr (E == EnergyType::ONE_SIDED_FIXED_JOINT)
            return one_sided_fixed_joint;
    }

    template<StaticEnergyType E>
    auto& get()
    {
        if constexpr (E == StaticEnergyType::NEO_HOOKEAN)
            return neo_hookean;
        else if constexpr (E == StaticEnergyType::COSSERAT_ROD)
            return cosserat_rod;
        else if constexpr (E == StaticEnergyType::GROUND_COLLISION)
            return ground_collision;
        else if constexpr (E == StaticEnergyType::RIGID_BODY_GROUND_COLLISION)
            return rigid_body_ground_collision;
        else if constexpr (E == StaticEnergyType::ONE_SIDED_FIXED_JOINT)
            return one_sided_fixed_joint;
    }

    template<DynamicEnergyType E>
    auto& get()
    {
        if constexpr (E == DynamicEnergyType::TRIANGLE_RIGID_COLLISION)
            return triangle_rigid_collision;
        else if constexpr (E == DynamicEnergyType::TRIANGLE_ROD_COLLISION)
            return triangle_rod_collision;
    }

    /** Apply a function to each set of energies. */
    template <typename Func>
    void forEachEnergyType(Func&& f)
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
        f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }
    template <typename Func>
    void forEachEnergyType(Func&& f) const
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
        f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }

    /** Apply a function only to energies derived from a hard constraint */
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f)
    {
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
        f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f) const
    {
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
        f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }

    /** Apply a function only to "static" energies which are generally constant throughout the simulation.
     * Examples: constitutive energies (NeoHookean, Corotational, etc.), Cosserat rod energies, etc.
     * 
     * This distinction between "static" and "dynamic" is ultimately which energies are handled by the seldom-rebuilt ParticleAdjacency
     * struct vs. those that are handled by the CollisionParticleAdjacency struct (which is rebuilt every frame)
     */
    template <typename Func>
    void forEachStaticEnergyType(Func&& f) 
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(one_sided_fixed_joint);
    }
    template <typename Func>
    void forEachStaticEnergyType(Func&& f) const
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(one_sided_fixed_joint);
    }

    /** Apply a function only to "dynamic" energies which are transient during the simulation.
     * Examples: dynamically generated collision constraints (triangle-triangle, triangle-rigid, etc.)
     * 
     * This distinction between "static" and "dynamic" is ultimately which energies are handled by the seldom-rebuilt ParticleAdjacency
     * struct vs. those that are handled by the CollisionParticleAdjacency struct (which is rebuilt every frame)
     */
    template <typename Func>
    void forEachDynamicEnergyType(Func&& f)
    {
        f(triangle_rigid_collision);
        f(triangle_rod_collision);
    }
    template <typename Func>
    void forEachDynamicEnergyType(Func&& f) const
    {
        f(triangle_rigid_collision);
        f(triangle_rod_collision);
    }
};

} // namespace Energy