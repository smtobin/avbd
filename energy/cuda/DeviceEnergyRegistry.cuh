#pragma once

#include "energy/EnergyRegistry.hpp"
#include "energy/cuda/DeviceNeoHookeanEnergyPool.cuh"
#include "energy/cuda/DeviceCosseratRodEnergyPool.cuh"
#include "energy/cuda/DeviceGroundCollisionEnergyPool.cuh"
#include "energy/cuda/DeviceRigidBodyGroundCollisionEnergyPool.cuh"
#include "energy/cuda/DeviceOneSidedFixedJointEnergyPool.cuh"

namespace Energy
{

struct DeviceEnergyRegistry
{
    // elastic
    DeviceNeoHookeanEnergyPool neo_hookean;
    DeviceCosseratRodEnergyPool cosserat_rod;

    // collision
    DeviceGroundCollisionEnergyPool ground_collision;
    DeviceRigidBodyGroundCollisionEnergyPool rigid_body_ground_collision;
    /** TODO: (10/23/26) Dynamic energies
     * 
     * Make sure to update member functions as well!
     */
    // DeviceTriangleRigidCollisionEnergyPool
    // DeviceTriangleRodCollisionEnergyPool

    // joints
    DeviceOneSidedFixedJointEnergyPool one_sided_fixed_joint;

    DeviceEnergyRegistry(unsigned capacity)
        : neo_hookean(capacity)
        , cosserat_rod(capacity)
        , ground_collision(capacity)
        , rigid_body_ground_collision(capacity)
        , one_sided_fixed_joint(capacity)
    {}

    void upload(const EnergyRegistry& src)
    {
        neo_hookean.upload(src.neo_hookean);
        cosserat_rod.upload(src.cosserat_rod);
        ground_collision.upload(src.ground_collision);
        rigid_body_ground_collision.upload(src.rigid_body_ground_collision);
        one_sided_fixed_joint.upload(src.one_sided_fixed_joint);
    }

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

    // template<DynamicEnergyType E>
    // auto& get()
    // {
    //     if constexpr (E == DynamicEnergyType::TRIANGLE_RIGID_COLLISION)
    //         return triangle_rigid_collision;
    //     else if constexpr (E == DynamicEnergyType::TRIANGLE_ROD_COLLISION)
    //         return triangle_rod_collision;
    // }

    /** Apply a function to each set of energies. */
    template <typename Func>
    void forEachEnergyType(Func&& f)
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        // f(triangle_rigid_collision);
        // f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }
    template <typename Func>
    void forEachEnergyType(Func&& f) const
    {
        f(neo_hookean);
        f(cosserat_rod);
        f(ground_collision);
        f(rigid_body_ground_collision);
        // f(triangle_rigid_collision);
        // f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }

    /** Apply a function only to energies derived from a hard constraint */
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f)
    {
        f(ground_collision);
        f(rigid_body_ground_collision);
        // f(triangle_rigid_collision);
        // f(triangle_rod_collision);
        f(one_sided_fixed_joint);
    }
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f) const
    {
        f(ground_collision);
        f(rigid_body_ground_collision);
        // f(triangle_rigid_collision);
        // f(triangle_rod_collision);
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
};

} // namespace Energy