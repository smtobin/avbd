#pragma once

#include "common/common.hpp"
#include "simulation/cuda/DeviceSimulationContext.cuh"
#include "simulation/cuda/DeviceVBDSolver.hpp"

#include <thread>

namespace Sim
{

/** High level management class
 * Manages members that launch kernels that performs all simulation work.
 * During each step:
 *  === Collision ===
 *  - construct BVH                        
 *  - traverse BVH for potential collisions
 *  - narrow-phase on potential collisions 
 *  - rebuild adjacency                    
 *  - rebuild color list                   
 *  === Solver Step ===
 *  - particle inertial update
 *  - initialize constraints
 *    === Solver Iteration ===
 *       - solve colors
 *       - Chebyshev acceleration
 *       - update constraints
 *  - velocity update
 */
class DeviceSimulationExecutor
{

private:
    /** Simulation context */
    DeviceSimulationContext* _ctx;

    /** Responsible for detecting collisions 
     * 
     * TODO: (10/23/26) Add this to SimulationExecutor, when its ready
    */
    // Collision::CollisionDetector _collision_detector;

    /** Responsible for implementing the AVBD algorithm */
    DeviceVBDSolver _solver;

    /** Number of threads */
    unsigned _num_threads;

public:
    DeviceSimulationExecutor() = default;

    /** Constructor - set up the worker threads */
    DeviceSimulationExecutor(
        DeviceSimulationContext* ctx,
        unsigned solver_iters,
        Real iter_acceleration,
        unsigned num_threads
    )
        : _ctx(ctx)
        // , _collision_detector(1000)
        , _solver(ctx, solver_iters, iter_acceleration)
        , _num_threads(num_threads)
    {
    }

    /** Destructor */
    ~DeviceSimulationExecutor() = default;

    void timeStep()
    {
        
    }

};

} // namespace Sim