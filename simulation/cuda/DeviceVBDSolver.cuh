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
    DeviceVBDSolver(DeviceSimulationContext* ctx, unsigned solver_iters, Real iter_acceleration);
    ~DeviceVBDSolver() = default;

    /** Launches kernels for a single time step */
    /** TODO: (10/23/26) write solver kernels */
    void solve();
};

} // namespace Sim