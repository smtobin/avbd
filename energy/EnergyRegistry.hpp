#pragma once

#include "energy/NeoHookeanEnergyPool.hpp"
#include "energy/GroundCollisionEnergyPool.hpp"
#include "energy/TriangleRigidCollisionEnergyPool.hpp"

namespace Energy
{

/** Storage of all the energies in the sim.
 * Stores the individual energy memory pools.
 */
struct EnergyRegistry
{
    NeoHookeanEnergyPool neo_hookean;
    GroundCollisionEnergyPool ground_collision;
    TriangleRigidCollisionEnergyPool triangle_rigid_collision;

    EnergyRegistry(unsigned capacity)
        : neo_hookean(capacity)
        , ground_collision(capacity)
        , triangle_rigid_collision(capacity)
    {}

    /** Statically-typed getter for a specific energy type */
    template<EnergyType E>
    auto& get()
    {
        if constexpr (E == EnergyType::NEO_HOOKEAN)
            return neo_hookean;
        else if constexpr (E == EnergyType::GROUND_COLLISION)
            return ground_collision;
        else if constexpr (E == EnergyType::TRIANGLE_RIGID_COLLISION)
            return triangle_rigid_collision;
    }

    /** Apply a function to each set of energies. */
    template <typename Func>
    void forEachEnergyType(Func&& f)
    {
        f(neo_hookean);
        f(ground_collision);
        // f(triangle_rigid_collision);
    }
    template <typename Func>
    void forEachEnergyType(Func&& f) const
    {
        f(neo_hookean);
        f(ground_collision);
        // f(triangle_rigid_collision);
    }

    /** Apply a function only to energies derived from a hard constraint */
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f)
    {
        f(ground_collision);
        f(triangle_rigid_collision);
    }
    template <typename Func>
    void forEachHardConstraintEnergyType(Func&& f) const
    {
        f(ground_collision);
        f(triangle_rigid_collision);
    }
};

} // namespace Energy