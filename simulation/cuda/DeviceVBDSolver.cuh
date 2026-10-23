#pragma once

#include "simulation/cuda/DeviceSimulationContext.cuh"

namespace Sim
{

class DeviceVBDSolver
{
private:
    /** Device simulation context */
    DeviceSimulationContext* _ctx;

    /** Number of solver iterations per time step */
    unsigned _solver_iters;

    /** Chebyshev acceleration parameter in [0, 1] */
    Real _iter_acceleration;

public:
    DeviceVBDSolver() = default;
    DeviceVBDSolver(DeviceSimulationContext* ctx, unsigned solver_iters, Real iter_acceleration)
        : _ctx(ctx)
        , _solver_iters(solver_iters)
        , _iter_accleration(iter_acceleration)
    {
    }

    ~DeviceVBDSolver() = default;

    /** Launches kernels for a single time step */
    /** TODO: (10/23/26) write solver kernels */
    void solve()
    {
        // update all energies after time step
        _ctx->energies.forEachEnergyType([&] (auto& pool) {
            // updateAfterTimeStep
        });
        
        
        // compute inertial update and initialization for each particle
        // particleInertialUpdate

        // solve each individual vertex block
        Real omega = 1;
        for (unsigned i = 0; i < _solver_iters; i++)
        {
            if (i > 2)
                omega = 4 / (4 - _iter_acceleration*_iter_acceleration*omega);
            else if (i == 2)
                omega = 2 / (2 - _iter_acceleration*_iter_acceleration);
            else
                omega = 1;

            // iterate through colors and parallelize within the color
            /** TODO: (10/23/26) num_colors currently lives on the device. Move it to host?  */
            for (unsigned c = 0; c < _ctx->coloring.num_colors; c++)
            {
                // solveParticleRangeinColor
            }

            // Chebyshev acceleration
            // particleChebyshevAcceleration

            // update Lagrange multipliers
            _ctx->energies.forEachHardConstraintEnergyType([&] (auto& pool) {
                // updateAfterIteration
            });
        }

        // velocity update
        // particleVelocityUpdate
    }
};

} // namespace Sim