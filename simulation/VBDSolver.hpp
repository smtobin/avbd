#pragma once

#include "common/common.hpp"
#include "common/SpinBarrier.hpp"
#include "common/EnergyUtils.hpp"
#include "simulation/SimulationContext.hpp"

#include "energy/NeoHookeanEnergySolver.hpp"
#include "energy/GroundCollisionEnergySolver.hpp"

#include <chrono>

namespace Sim
{

class VBDSolver
{
private:
    /** Simulation context */
    SimulationContext* _ctx;  

    /** Number of solver iterations per time step */
    unsigned _solver_iters;

    /** Chebyshev acceleration parameter in [0, 1] */
    Real _iter_acceleration;

    /** Number of threads */
    unsigned _num_threads;

    /** Worker threads */
    std::vector<std::thread> _workers;

    /** If the workers are actively running */
    std::atomic<bool> _running;
    
    /** Lightweight synchronization barrier for each worker thread */
    SpinBarrier _barrier;

    /** Solver generation used as start signal for worker threads */
    std::atomic<unsigned> _solver_generation;

public:
    VBDSolver() = default;

    VBDSolver(SimulationContext* ctx, unsigned solver_iters, Real iter_acceleration, unsigned num_threads)
        : _ctx(ctx)
        , _solver_iters(solver_iters)
        , _iter_acceleration(iter_acceleration)
        , _num_threads(num_threads)
        , _running(true)
        , _barrier(num_threads)
        , _solver_generation(0)
    {
        // create the worker threads
        _workers.reserve(num_threads);
        for (unsigned w_idx = 1; w_idx < num_threads; w_idx++)
        {
            _workers.emplace_back([this, w_idx] {_workerThread(w_idx); });
        }
    }

    ~VBDSolver()
    {
        _running.store(false);

        _solver_generation.fetch_add(1, std::memory_order_release);

        for (auto& t : _workers)
            t.join();
    }

    void _workerThread(unsigned w_idx)
    {
        unsigned this_gen = _solver_generation.load();
        while (_running)
        {
            // wait for start signal
            unsigned test_gen;
            do
            {
                test_gen = _solver_generation.load(std::memory_order_acquire);
                _mm_pause();
            } while (test_gen == this_gen);
            this_gen = test_gen;

            _workerIteration(w_idx);
        }
    }

    void _workerIteration(unsigned w_idx)
    {   
        Vec3r a_grav(0, -_ctx->params.g_accel, 0);
        Real dt = _ctx->params.dt;

        // compute inertial update and initialization for each particle
        _particleRangeInertialUpdate(w_idx, a_grav, dt);
        _barrier.arrive_and_wait();

        // update all energies after time step
        _ctx->energies.forEachEnergyType([&] (auto& pool) {
            using Pool = base_type_t<decltype(pool)>;

            _updateRangeAfterTimeStep(w_idx, pool, _ctx->particles);
            _barrier.arrive_and_wait();
        });

        // solve each individual vertex block
        Real omega = 1;
        for (unsigned i = 0; i < _solver_iters; i++)
        {
            // std::cout << "=== Starting iter " << i << " === " << std::endl;
            if (i > 2)
                omega = 4 / (4 - _iter_acceleration*_iter_acceleration*omega);
            else if (i == 2)
                omega = 2 / (2 - _iter_acceleration*_iter_acceleration);
            else
                omega = 1;

            // iterate through colors and parallelize within the color
            for (unsigned c = 0; c < _ctx->coloring.num_colors; c++)
            {   
                _solveParticleRangeInColor(w_idx, c, dt);
                _barrier.arrive_and_wait();

                // copy buffer into vertices
                // _ctx->particles.positions = _ctx->particles.buffered_positions;
            }

            // Chebyshev acceleration
            _particleRangeChebyshevAcceleration(w_idx, omega);
            _barrier.arrive_and_wait();

            // update Lagrange multipliers
            _ctx->energies.forEachHardConstraintEnergyType([&] (auto& pool) {
                using Pool = base_type_t<decltype(pool)>;

                _updateRangeAfterIteration(w_idx, pool, _ctx->particles);
                _barrier.arrive_and_wait();
            });
        }

        _particleRangeVelocityUpdate(w_idx, dt);
        _barrier.arrive_and_wait();
    }

    /** Compute start and end indices for worker thread when iterating over a TombstonePool. */
    std::pair<unsigned, unsigned> _computeStartEnd(unsigned w_idx, const TombstonePool& pool)
    {
        unsigned num = pool.highest_index+1;
        return _computeStartEnd(w_idx, pool.highest_index+1);
    }

    /** Compute start and end indices for worker thread when iterating over a number of objects. */
    std::pair<unsigned, unsigned> _computeStartEnd(unsigned w_idx, unsigned total_num)
    {
        unsigned start = total_num * w_idx / _num_threads;
        unsigned end   = total_num * (w_idx + 1) / _num_threads;
        return std::make_pair(start, end);
    }

    /** Worker thread inertial update over its range of particles */
    void _particleRangeInertialUpdate(unsigned w_idx, const Vec3r& a_ext, Real dt)
    {
        auto [start, end] = _computeStartEnd(w_idx, _ctx->particles);
        for (unsigned p_idx = start; p_idx < end; p_idx++)
        {
            _particleInertialUpdate(p_idx, a_ext, dt);
        }
    }

    template <typename EnergyPool>
    void _updateRangeAfterTimeStep(unsigned w_idx, EnergyPool& pool, ParticlePool& particles)
    {
        auto [start, end] = _computeStartEnd(w_idx, pool);
        for (unsigned e_idx = start; e_idx < end; e_idx++)
        {
            if (pool.active[e_idx])
                EnergyPool::Solver::updateAfterTimeStep(e_idx, pool, particles);
        }
    }

    template <typename EnergyPool>
    void _updateRangeAfterIteration(unsigned w_idx, EnergyPool& pool, ParticlePool& particles)
    {
        auto [start, end] = _computeStartEnd(w_idx, pool);
        for (unsigned e_idx = start; e_idx < end; e_idx++)
        {
            if (pool.active[e_idx])
                EnergyPool::Solver::updateAfterIteration(e_idx, pool, particles);
        }
    }

    void _solveParticleRangeInColor(unsigned w_idx, unsigned color, Real dt)
    {
        auto [start, end] = _computeStartEnd(w_idx, _ctx->coloring.color_counts[color]);
        unsigned color_start = _ctx->coloring.color_offsets[color];
        for (unsigned c_idx = color_start + start; c_idx < color_start + end; c_idx++)
        {
            unsigned p_idx = _ctx->coloring.work_list[c_idx];
            _solveParticle(p_idx, dt);
        }
    }

    void _particleRangeChebyshevAcceleration(unsigned w_idx, Real omega)
    {
        auto [start, end] = _computeStartEnd(w_idx, _ctx->particles);
        for (unsigned p_idx = start; p_idx < end; p_idx++)
        {
            if (_ctx->particles.active[p_idx])
                _particleChebyshevAcceleration(p_idx, omega);
        }
    }

    void _particleRangeVelocityUpdate(unsigned w_idx, Real dt)
    {
        auto [start, end] = _computeStartEnd(w_idx, _ctx->particles);
        for (unsigned p_idx = start; p_idx < end; p_idx++)
        {
            if (_ctx->particles.active[p_idx])
                _particleVelocityUpdate(p_idx, dt);
        }
    }

    void solve(Real dt)
    {
        // increment solver generation to signal worker threads to start
        _solver_generation.fetch_add(1, std::memory_order_release);

        // use this thread to perform a worker iteration (this is thread 0)
        // synchronization at the end of the worker iteration guarantees that we finish the work before moving on
        _workerIteration(0);
    }

private:
    void _particleInertialUpdate(unsigned p_idx, const Vec3r& a_ext, Real dt)
    {
        // std::cout << "\n=== Particle " << p_idx << " inertial update" << std::endl;
        Vec3r& p = _ctx->particles.positions[p_idx];
        Vec3r& y = _ctx->particles.inertial_positions[p_idx];
        const Vec3r& v = _ctx->particles.velocities[p_idx];
        const Vec3r& v_prev = _ctx->particles.previous_velocities[p_idx];

        // use adaptive initialization (Sec 3.7 in VBD paper)
        Vec3r a = (v - v_prev) / dt;

        Vec3r a_tilde_vec = Vec3r::Zero();
        Real a_ext_norm = a_ext.norm();
        if (a_ext_norm > 1e-12)
        {
            Real a_along_a_ext = a.dot(a_ext) / a_ext_norm;

            Real a_tilde;
            if (a_along_a_ext > a_ext_norm)
                a_tilde = 1;
            else if (a_along_a_ext < 0)
                a_tilde = 0;
            else
                a_tilde = a_along_a_ext / a_ext_norm;

            a_tilde_vec = a_tilde * a_ext;
        }

        // compute the inertially predicted position
        y = p + dt*v + dt*dt*a_ext;

        // move particle to its initialized position
        p += dt*v + dt*dt*a_tilde_vec;

        // mark particle as not in collision at the beginning of the time step
        _ctx->particles.in_collision[p_idx] = false;

        // initialize the previous iteration positions
        _ctx->particles.last_iter_positions[p_idx] = _ctx->particles.positions[p_idx];
        _ctx->particles.last_last_iter_positions[p_idx] = _ctx->particles.positions[p_idx];
    }

    void _solveParticle(unsigned p_idx, Real dt)
    {
        // std::cout << "\n=== Particle " << p_idx << " solve" << std::endl;

        const PerEnergy<unsigned>& adj_offsets = _ctx->adjacency.e_offsets[p_idx];

        Vec3r grad = Vec3r::Zero();
        Mat3r hess = Mat3r::Zero();

        // iterate through energy types
        unsigned num_energies = static_cast<unsigned>(EnergyType::size);
        ForEachEnergy([&]<EnergyType E>() {
            using Solver = SolverFor<E>::type;
            if constexpr (Energy::HasAccumulate4<Solver>)
            {
                // accumulate 4 loop
                // Solver::accumulate4(...);
            }
            else
            {
                // normal single loop
            }
        }
        );

        
        

        // std::cout << "=== Starting Accumulate ===" << std::endl;
        /**  TODO: this is hard-coded for the case of N NeoHookean energies and 1 GroundCollision, which is the case for the sim currently.
         * Need to generalize this to more constraints that are not necessarily arranged in this structure.
         */
        unsigned e = adj_offsets[0];
        unsigned adj_end = _ctx->adjacency.e_offsets[p_idx+1][0];    // the end of the range is the start of the next offsets
        for (; e+3 < adj_end-1; e+=4)
        {
            const ParticleAdjacency::Entry& entry1 = _ctx->adjacency.e_entries[e];
            const ParticleAdjacency::Entry& entry2 = _ctx->adjacency.e_entries[e+1];
            const ParticleAdjacency::Entry& entry3 = _ctx->adjacency.e_entries[e+2];
            const ParticleAdjacency::Entry& entry4 = _ctx->adjacency.e_entries[e+3];
            unsigned e_idx[4] = {entry1.energy_idx,
                entry2.energy_idx,
                entry3.energy_idx,
                entry4.energy_idx};
            unsigned l_idx[4] = {entry1.local_vertex_idx,
                entry2.local_vertex_idx,
                entry3.local_vertex_idx,
                entry4.local_vertex_idx};
            Energy::NeoHookeanEnergySolver::accumulate4(
                e_idx,
                _ctx->energies.neo_hookean,
                _ctx->particles,
                l_idx,
                hess,
                grad,
                dt
            );
        }
        for (; e < adj_end-1; e++)
        {
            const ParticleAdjacency::Entry& entry = _ctx->adjacency.e_entries[e];
            Energy::NeoHookeanEnergySolver::accumulate(
                entry.energy_idx,
                _ctx->energies.neo_hookean,
                _ctx->particles,
                entry.local_vertex_idx,
                hess,
                grad,
                dt
            );
        }
        const ParticleAdjacency::Entry& entry = _ctx->adjacency.e_entries[e];
        if (entry.energy_type == EnergyType::GROUND_COLLISION) 
        {
            // std::cout << " GroundCollision constraint " << entry.energy_idx << std::endl;
            Energy::GroundCollisionEnergySolver::accumulate(
                entry.energy_idx,
                _ctx->energies.ground_collision,
                _ctx->particles,
                entry.local_vertex_idx,
                hess,//H_acc[e & 3],
                grad,//G_acc[e & 3],
                dt
            );
        }

        // assemble LHS and RHS of single-particle system
        Real mass = _ctx->particles.masses[p_idx];
        const Vec3r& p = _ctx->particles.positions[p_idx];
        const Vec3r& y = _ctx->particles.inertial_positions[p_idx];

        Vec3r RHS = -mass / (dt*dt) * (p - y) - grad;
        Mat3r LHS = mass / (dt*dt) * Mat3r::Identity() + hess;

        // Vec3r dx = LHS.partialPivLu().solve(RHS);
        Vec3r dx = LHS.inverse() * RHS;

        // std::cout << "dx: " << dx.transpose() << std::endl;

        _ctx->particles.positions[p_idx] += dx;
        // put positions into a buffer
        // _ctx->particles.buffered_positions[p_idx] = _ctx->particles.positions[p_idx] + dx;
    }

    void _particleChebyshevAcceleration(unsigned p_idx, Real omega)
    {
        // only do Chebyshev acceleration for particles not in collision
        if (!_ctx->particles.in_collision[p_idx])
        {
            _ctx->particles.positions[p_idx] =
                omega * (_ctx->particles.positions[p_idx] - _ctx->particles.last_last_iter_positions[p_idx]) + _ctx->particles.last_last_iter_positions[p_idx];
        }

        // update the previous iteration positions
        _ctx->particles.last_last_iter_positions[p_idx] = _ctx->particles.last_iter_positions[p_idx];
        _ctx->particles.last_iter_positions[p_idx] = _ctx->particles.positions[p_idx];
    }

    void _particleVelocityUpdate(unsigned p_idx, Real dt)
    {
        _ctx->particles.previous_velocities[p_idx] = _ctx->particles.velocities[p_idx];
        _ctx->particles.velocities[p_idx] = 
            (_ctx->particles.positions[p_idx] - _ctx->particles.previous_positions[p_idx]) / dt;

        _ctx->particles.previous_positions[p_idx] = _ctx->particles.positions[p_idx];
    }
};

} // namespace Sim