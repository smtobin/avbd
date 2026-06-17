#pragma once

#include "common/common.hpp"
#include "simulation/SimulationContext.hpp"

#include "energy/NeoHookeanEnergySolver.hpp"
#include "energy/GroundCollisionEnergySolver.hpp"

namespace Sim
{

class VBDSolver
{
private:
    /** Simulation context */
    SimulationContext* _ctx;  

public:
    VBDSolver(SimulationContext* ctx)
        : _ctx(ctx)
    {

    }

    void solve(Real dt)
    {
        // iterate through particles and solve system
        /** TODO: When graph coloring is used, we will iterate through all the energies in the color. */
        for (unsigned p_idx : _ctx->particles)
        {
            _solveParticle(p_idx, dt);
        }
    }

private:
    void _solveParticle(unsigned p_idx, Real dt)
    {
        unsigned adj_start = _ctx->adjacency.adjOffsets[p_idx];
        unsigned adj_end   = _ctx->adjacency.adjOffsets[p_idx + 1];
        for (unsigned e = adj_start; e < adj_end; e++) 
        {
            const ParticleAdjacency::Entry& entry = _ctx->adjacency.adj_entries[e];
            if (entry.energy_type == EnergyType::NEO_HOOKEAN)
            {
                // Energy::NeoHookeanEnergySolver::accumulate()
            } 
            else if (entry.energy_type == EnergyType::GROUND_COLLISION) 
            {
                // Energy::GroundCollisionEnergySolver::accumulate()
            }
        }
    }
};

} // namespace Sim