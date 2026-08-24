#pragma once

#include "energy/NeoHookeanEnergyPool.hpp"
#include "energy/CosseratRodEnergyPool.hpp"
#include "energy/GroundCollisionEnergyPool.hpp"
#include "energy/RigidBodyGroundCollisionEnergyPool.hpp"
#include "energy/TriangleRigidCollisionEnergyPool.hpp"

namespace Energy
{

/** Storage of all the energies in the sim.
 * Stores the individual energy memory pools.
 */
struct EnergyRegistry
{
    NeoHookeanEnergyPool neo_hookean;
    CosseratRodEnergyPool cosserat_rod;
    GroundCollisionEnergyPool ground_collision;
    RigidBodyGroundCollisionEnergyPool rigid_body_ground_collision;
    TriangleRigidCollisionEnergyPool triangle_rigid_collision;

    EnergyRegistry(unsigned capacity)
        : neo_hookean(capacity)
        , cosserat_rod(capacity)
        , ground_collision(capacity)
        , rigid_body_ground_collision(capacity)
        , triangle_rigid_collision(capacity)
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
    }

    template<DynamicEnergyType E>
    auto& get()
    {
        if constexpr (E == DynamicEnergyType::TRIANGLE_RIGID_COLLISION)
            return triangle_rigid_collision;
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
    }
    template <typename Func>
    void forEachEnergyType(Func&& f) const
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
    }

    /** Apply a function only to energies derived from a hard constraint */
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f)
    {
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
    }
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f) const
    {
        f(ground_collision);
        f(rigid_body_ground_collision);
        f(triangle_rigid_collision);
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
    }
    template <typename Func>
    void forEachStaticEnergyType(Func&& f) const
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
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
    }
    template <typename Func>
    void forEachDynamicEnergyType(Func&& f) const
    {
        f(triangle_rigid_collision);
    }
};

} // namespace Energy