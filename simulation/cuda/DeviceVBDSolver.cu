#include "simulation/cuda/DeviceVBDSolver.cuh"

namespace Sim
{

namespace
{

// Kernel definitions private to this implementation file.

__global__ void particleInertialUpdateKernel(/* device data arguments */)
{
    // One thread per particle, for example.
}

__global__ void particleChebyshevAccelerationKernel(/* device data arguments */)
{
    // Update particle positions using the acceleration parameter.
}

__global__ void particleVelocityUpdateKernel(/* device data arguments */)
{
    // Update particle velocities.
}

} // anonymous namespace

DeviceVBDSolver::DeviceVBDSolver(
    DeviceSimulationContext* ctx,
    unsigned solver_iters,
    Real iter_acceleration)
    : _ctx(ctx)
    , _solver_iters(solver_iters)
    , _iter_acceleration(iter_acceleration)
{
}

void DeviceVBDSolver::solve()
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

} // namespace Sim