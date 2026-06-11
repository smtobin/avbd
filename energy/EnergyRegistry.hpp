#pragma once

#include "energy/NeoHookeanEnergyPool.hpp"
#include "energy/GroundCollisionEnergyPool.hpp"

namespace Energy
{

/** Storage of all the energies in the sim.
 * Stores the individual energy memory pools.
 */
struct EnergyRegistry
{
    NeoHookeanEnergyPool neo_hookean;
    GroundCollisionEnergyPool ground_collision;

    /** Apply a function to each set of energies. */
    template <typename Func>
    void forEachEnergyType(Func&& f)
    {
        f(neo_hookean);
        f(ground_collision);
    }
};

} // namespace Energy