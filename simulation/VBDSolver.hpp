#pragma once

#include "common/common.hpp"
#include "common/SpinBarrier.hpp"
#include "common/EnergyUtils.hpp"
#include "simulation/SimulationContext.hpp"

#include "energy/NeoHookeanEnergySolver.hpp"
#include "energy/GroundCollisionEnergySolver.hpp"
#include "energy/TriangleRigidCollisionEnergySolver.hpp"

#include <chrono>
#include <thread>

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

public:
    VBDSolver() = default;

    VBDSolver(SimulationContext* ctx, unsigned solver_iters, Real iter_acceleration)
        : _ctx(ctx)
        , _solver_iters(solver_iters)
        , _iter_acceleration(iter_acceleration)
    {
    }

    ~VBDSolver()
    {
    }

    void solve_Parallel(WorkerThreadContext& w_ctx)
    {   
        Vec3r a_grav(0, -_ctx->params.g_accel, 0);
        Real dt = _ctx->params.dt;

        // compute inertial update and initialization for each particle
        _particleRangeInertialUpdate(w_ctx, a_grav, dt);
        w_ctx.barrier->arrive_and_wait();

        // update all energies after time step
        _ctx->energies.forEachEnergyType([&] (auto& pool) {
            using Pool = base_type_t<decltype(pool)>;

            _updateRangeAfterTimeStep(w_ctx, pool, _ctx->particles);
            w_ctx.barrier->arrive_and_wait();
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
                _solveParticleRangeInColor(w_ctx, c, dt);
                w_ctx.barrier->arrive_and_wait();

                // only pay the extra sync when this color actually has leftover conflicts
                if (_ctx->coloring.conflicted_by_color_offsets[c+1] != _ctx->coloring.conflicted_by_color_offsets[c])
                {
                    _commitConflictedInColor(w_ctx, c);
                    w_ctx.barrier->arrive_and_wait();
                }
            }

            // Chebyshev acceleration
            _particleRangeChebyshevAcceleration(w_ctx, omega);
            w_ctx.barrier->arrive_and_wait();

            // update Lagrange multipliers
            _ctx->energies.forEachHardConstraintEnergyType([&] (auto& pool) {
                using Pool = base_type_t<decltype(pool)>;

                _updateRangeAfterIteration(w_ctx, pool, _ctx->particles);
                w_ctx.barrier->arrive_and_wait();
            });
        }

        _particleRangeVelocityUpdate(w_ctx, dt);
        w_ctx.barrier->arrive_and_wait();
    }

    /** Commits particles within a color that are conflicted.
     * During the solve step, conflicted particle updates are buffered, so we must commit the buffered positions.
     */
    void _commitConflictedInColor(WorkerThreadContext& w_ctx, unsigned c)
    {
        unsigned start = _ctx->coloring.conflicted_by_color_offsets[c];
        unsigned end   = _ctx->coloring.conflicted_by_color_offsets[c+1];

        auto [lo, hi] = w_ctx.computeStartEnd(end - start);
        for (unsigned k = start + lo; k < start + hi; k++)
        {
            unsigned p = _ctx->coloring.conflicted_by_color_entries[k];
            _ctx->particles.positions[p] = _ctx->particles.buffered_positions[p];
        }
    }

    /** Worker thread inertial update over its range of particles */
    void _particleRangeInertialUpdate(WorkerThreadContext& w_ctx, const Vec3r& a_ext, Real dt)
    {
        auto [start, end] = w_ctx.computeStartEnd(_ctx->particles);
        for (unsigned p_idx = start; p_idx < end; p_idx++)
        {
            _particleInertialUpdate(p_idx, a_ext, dt);
        }
    }

    template <typename EnergyPool>
    void _updateRangeAfterTimeStep(WorkerThreadContext& w_ctx, EnergyPool& pool, ParticlePool& particles)
    {
        auto [start, end] = w_ctx.computeStartEnd(pool);
        for (unsigned e_idx = start; e_idx < end; e_idx++)
        {
            if (pool.active[e_idx])
                EnergyPool::SolverType::updateAfterTimeStep(e_idx, pool, particles);
        }
    }

    template <typename EnergyPool>
    void _updateRangeAfterIteration(WorkerThreadContext& w_ctx, EnergyPool& pool, ParticlePool& particles)
    {
        auto [start, end] = w_ctx.computeStartEnd(pool);
        for (unsigned e_idx = start; e_idx < end; e_idx++)
        {
            if (pool.active[e_idx])
                EnergyPool::SolverType::updateAfterIteration(e_idx, pool, particles);
        }
    }

    void _solveParticleRangeInColor(WorkerThreadContext& w_ctx, unsigned color, Real dt)
    {
        auto [start, end] = w_ctx.computeStartEnd(_ctx->coloring.color_counts[color]);
        unsigned color_start = _ctx->coloring.color_offsets[color];
        for (unsigned c_idx = color_start + start; c_idx < color_start + end; c_idx++)
        {
            unsigned p_idx = _ctx->coloring.work_list[c_idx];
            _solveParticle(p_idx, dt);
        }
    }

    void _particleRangeChebyshevAcceleration(WorkerThreadContext& w_ctx, Real omega)
    {
        auto [start, end] = w_ctx.computeStartEnd(_ctx->particles);
        for (unsigned p_idx = start; p_idx < end; p_idx++)
        {
            if (_ctx->particles.active[p_idx])
                _particleChebyshevAcceleration(p_idx, omega);
        }
    }

    void _particleRangeVelocityUpdate(WorkerThreadContext& w_ctx, Real dt)
    {
        auto [start, end] = w_ctx.computeStartEnd(_ctx->particles);
        for (unsigned p_idx = start; p_idx < end; p_idx++)
        {
            if (_ctx->particles.active[p_idx])
                _particleVelocityUpdate(p_idx, dt);
        }
    }

    // void solve(Real dt)
    // {
    //     // increment solver generation to signal worker threads to start
    //     _solver_generation.fetch_add(1, std::memory_order_release);

    //     // use this thread to perform a worker iteration (this is thread 0)
    //     // synchronization at the end of the worker iteration guarantees that we finish the work before moving on
    //     _workerIteration(0);
    // }

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
        if (_ctx->particles.isOriented(p_idx))
            return;
        // std::cout << "\n=== Particle " << p_idx << " solve" << std::endl;

        Vec3r grad = Vec3r::Zero();
        Mat3r hess = Mat3r::Zero();

        /** Process "static" energies */
        {
            const PerStaticEnergy<unsigned>& static_adj_offsets = _ctx->static_adjacency.e_offsets[p_idx];
            unsigned adj_end = _ctx->static_adjacency.e_offsets[p_idx+1][0];

            // iterate through energy types
            unsigned num_energies = static_cast<unsigned>(StaticEnergyType::count);
            Energy::ForEachStaticEnergy([&]<StaticEnergyType E>() {
                unsigned e_type_idx = static_cast<unsigned>(E);
                // get the starting and ending offsets from the adjacency list
                // the start index of the adjacent energies of this type
                unsigned e_adj = static_adj_offsets[e_type_idx]; 
                // the end index (use conditional to wrap around if needed)
                unsigned final_e_adj = (e_type_idx+1 == num_energies) ? adj_end : static_adj_offsets[e_type_idx+1];

                using Solver = SolverFor<E>::type;
                const auto& energy_pool = _ctx->energies.template get<E>();

                // if solver has AVX accumulate4 implemented, use that
                if constexpr (Energy::HasAccumulate4<Solver>)
                {
                    // process in chunks of 4 using AVX
                    for (; e_adj+3 < final_e_adj; e_adj+=4)
                    {
                        const StaticParticleAdjacency::Entry& entry1 = _ctx->static_adjacency.e_entries[e_adj];
                        const StaticParticleAdjacency::Entry& entry2 = _ctx->static_adjacency.e_entries[e_adj+1];
                        const StaticParticleAdjacency::Entry& entry3 = _ctx->static_adjacency.e_entries[e_adj+2];
                        const StaticParticleAdjacency::Entry& entry4 = _ctx->static_adjacency.e_entries[e_adj+3];
                        unsigned e_idx[4] = {entry1.energy_idx,
                            entry2.energy_idx,
                            entry3.energy_idx,
                            entry4.energy_idx};
                        unsigned l_idx[4] = {entry1.local_vertex_idx,
                            entry2.local_vertex_idx,
                            entry3.local_vertex_idx,
                            entry4.local_vertex_idx};
                        Solver::accumulate4(
                            e_idx,
                            energy_pool,
                            _ctx->particles,
                            l_idx,
                            hess,
                            grad,
                            dt
                        );
                    }

                    // process remainder
                    for (; e_adj < final_e_adj; e_adj++)
                    {
                        const StaticParticleAdjacency::Entry& entry = _ctx->static_adjacency.e_entries[e_adj];
                        Solver::accumulate(
                            entry.energy_idx,
                            energy_pool,
                            _ctx->particles,
                            entry.local_vertex_idx,
                            hess,
                            grad,
                            dt
                        );
                    }
                }
                // otherwise fall back to normal 1-by-1 computation
                else
                {
                    for (; e_adj < final_e_adj; e_adj++)
                    {
                        const StaticParticleAdjacency::Entry& entry = _ctx->static_adjacency.e_entries[e_adj];
                        Solver::accumulate(
                            entry.energy_idx,
                            energy_pool,
                            _ctx->particles,
                            entry.local_vertex_idx,
                            hess,
                            grad,
                            dt
                        );
                    }
                }
            }
            );
        }

        /** Process "dynamic" energies */
        {
            const PerDynamicEnergy<unsigned>& dyn_adj_offsets = _ctx->dynamic_adjacency.e_offsets[p_idx];
            unsigned adj_end = _ctx->dynamic_adjacency.e_offsets[p_idx+1][0];

            // iterate through energy types
            unsigned num_energies = static_cast<unsigned>(DynamicEnergyType::count);
            Energy::ForEachDynamicEnergy([&]<DynamicEnergyType E>() {
                unsigned e_type_idx = static_cast<unsigned>(E);
                // get the starting and ending offsets from the adjacency list
                // the start index of the adjacent energies of this type
                unsigned e_adj = dyn_adj_offsets[e_type_idx]; 
                // the end index (use conditional to wrap around if needed)
                unsigned final_e_adj = (e_type_idx+1 == num_energies) ? adj_end : dyn_adj_offsets[e_type_idx+1];

                using Solver = SolverFor<E>::type;
                const auto& energy_pool = _ctx->energies.template get<E>();
                // dynamic constraints always 1-by-1 computation
                for (; e_adj < final_e_adj; e_adj++)
                {
                    const DynamicParticleAdjacency::Entry& entry = _ctx->dynamic_adjacency.e_entries[e_adj];
                    Solver::accumulate(
                        entry.energy_idx,
                        energy_pool,
                        _ctx->particles,
                        entry.local_vertex_idx,
                        hess,
                        grad,
                        dt
                    );
                }
            }
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

        // if this particle has an intra-color conflict, buffer its update
        if (_ctx->coloring.is_conflicted[p_idx])
            _ctx->particles.buffered_positions[p_idx] = _ctx->particles.positions[p_idx] + dx;
        // otherwise handle it normally
        else
            _ctx->particles.positions[p_idx] += dx;
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