#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"
#include "energy/EnergyRegistry.hpp"

namespace Simulation
{

/** Storage of all state in the sim. */
struct SimulationContext
{
    // memory pools
    ParticlePool particles;
    Energy::EnergyRegistry energies;
};

} // namespace Simulation