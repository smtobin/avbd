#pragma once

#include "common/common.hpp"

namespace Energy
{

template<typename Solver>
concept HasAccumulate4 = requires(
    unsigned e_idx[4],
    const NeoHookeanEnergyPool& energies,
    ParticlePool& particles,
    unsigned local_idx[4],
    Mat3r& particle_H,
    Vec3r& particle_G,
    Real dt
)
{
    Solver::accumulate4(e_idx, energies, particles, local_idx, particle_H, particle_G, dt);
};

}