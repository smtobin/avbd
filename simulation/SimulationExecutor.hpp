#pragma once

#include "common/common.hpp"
#include "common/WorkerThreadContext.hpp"
#include "simulation/SimulationContext.hpp"
#include "simulation/VBDSolver.hpp"

#include <thread>

namespace Sim
{

/** Class that manages worker threads which ultimately performs all simulation work.
 * During each step, the threads:
 *  === Collision ===
 *  - construct BVH                         TODO: (07/30/26) make parallel
 *  - traverse BVH for potential collisions TODO: (07/30/26) make parallel
 *  - narrow-phase on potential collisions  TODO: (07/30/26) make parallel
 *  - rebuild adjacency                     TODO: (07/30/26) make parallel
 *  - rebuild color list                    TODO: (07/30/26) make parallel
 *  === Solver Step ===
 *  - particle inertial update
 *  - initialize constraints
 *    === Solver Iteration ===
 *       - solve colors
 *       - Chebyshev acceleration
 *       - update constraints
 *  - velocity update
 */
class SimulationExecutor
{

private:
    /** Simulation context */
    SimulationContext* _ctx;

    /** Responsible for detecting collisions */
    Collision::CollisionDetector _collision_detector;

    /** Responsible for implementing the AVBD algorithm */
    VBDSolver _solver;

    /** Number of threads */
    unsigned _num_threads;

    /** Worker threads */
    std::vector<std::thread> _workers;

    /** Memory for worker threads */
    std::vector<WorkerThreadContext> _worker_contexts;

    /** If the workers are actively running */
    std::atomic<bool> _running;
    
    /** Lightweight synchronization barrier for each worker thread */
    SpinBarrier _barrier;

    /** Solver generation used as start signal for worker threads */
    std::atomic<unsigned> _generation;

public:
    SimulationExecutor() = default;

    /** Constructor - set up the worker threads */
    SimulationExecutor(SimulationContext* ctx, unsigned solver_iters, Real iter_acceleration, unsigned num_threads)
        : _ctx(ctx)
        , _collision_detector(1000) /** TODO: (07/31/26) capacity for collision detector? */
        , _solver(ctx, solver_iters, iter_acceleration)
        , _num_threads(num_threads)
        , _running(true)
        , _barrier(num_threads)
        , _generation(0)
    {
        
        WorkerThreadContext::NUM_THREADS = num_threads;

        // create worker context for thread 0
        _worker_contexts.reserve(num_threads);
        _worker_contexts.emplace_back(0, &_barrier);

        // create worker threads
        _workers.reserve(num_threads);
        for (unsigned w_idx = 1; w_idx < num_threads; w_idx++)
        {
            _worker_contexts.emplace_back(w_idx, &_barrier);
            _workers.emplace_back([this, w_idx] {_workerThread(w_idx); });
        }
    }

    /** Delete copy and move since worker threads capture "this" */
    SimulationExecutor(const SimulationExecutor&) = delete;
    SimulationExecutor(SimulationExecutor&&) = delete;
    SimulationExecutor& operator= (const SimulationExecutor&) = delete;
    SimulationExecutor& operator= (SimulationExecutor&&) = delete;

    /** Destructor - shut down worker threads */
    ~SimulationExecutor()
    {
        _running.store(false);

        _generation.fetch_add(1, std::memory_order_release);

        for (auto& t : _workers)
            t.join();
    }

    void timeStep()
    {
        // increment generation to start the time step
        _generation.fetch_add(1, std::memory_order_release);

        // use this thread to perform a worker iteration (this is thread 0)
        // synchronization at the end of the worker iteration guarantees that we finish the work before moving on
        _workerTimeStep(0);
    }

    /** Each worker thread executes this function until the simulation ends.
     * The thread will spin, waiting for _generation to be incremented indicating that a new step should be taken.
     * _workerIteration() has barriers within it for thread-level synchronization.
     * 
     * @param w_idx : the worker thread index
     */
    void _workerThread(unsigned w_idx)
    {
        unsigned this_gen = _generation.load();
        while (_running)
        {
            // wait for start signal
            unsigned test_gen;
            do
            {
                test_gen = _generation.load(std::memory_order_acquire);
                _mm_pause();
            } while (test_gen == this_gen);
            this_gen = test_gen;

            _workerTimeStep(w_idx);
        }
    }

    /** A single time step in the sim, executed by each worker thread in parallel */
    void _workerTimeStep(unsigned w_idx)
    {
        /** Collision detection */
        _collision_detector.detectCollisionsAndRecolor_Parallel(_worker_contexts[w_idx], _worker_contexts, *_ctx);

        /** AVBD solver */
        _solver.solve_Parallel(_worker_contexts[w_idx]);
        
    }

};

} // namespace Sim