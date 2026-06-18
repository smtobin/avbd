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

    /** Number of solver iterations per time step */
    unsigned _solver_iters;

    /** Chebyshev acceleration parameter in [0, 1] */
    Real _iter_acceleration;

public:
    VBDSolver(SimulationContext* ctx, unsigned solver_iters, Real iter_acceleration)
        : _ctx(ctx)
        , _solver_iters(solver_iters)
        , _iter_acceleration(iter_acceleration)
    {

    }

    void solve(Real dt)
    {
        Vec3r a_grav(0, -9.81, 0);  // acceleration due to gravity - applied to all particles
        // compute inertial update and initialization for each particle
        for (unsigned p_idx : _ctx->particles)
        {
            _particleInertialUpdate(p_idx, a_grav, dt);
        }

        // for all hard constraint energies with Lagrange multipliers, reset the multipliers and stiffnesses for the new time step
        _ctx->energies.forEachHardConstraintEnergyType([&] (const auto& pool) {
            for (unsigned e_idx : pool)
            {
                /** TODO: this is hard-coded for now. Make type detectino automatic */
                Energy::HardConstraintEnergySolver<Energy::GroundCollisionEnergyPool, Energy::GroundCollisionConstraintSolver>::updateAfterTimeStep(
                    e_idx,
                    pool
                );
            }
        });

        // solve each individual vertex block
        Real omega = 1;
        for (int i = 0; i < _solver_iters; i++)
        {
            std::cout << "=== Starting iter " << i << " === " << std::endl;
            if (i > 2)
                omega = 4 / (4 - _iter_acceleration*_iter_acceleration*omega);
            else if (i == 2)
                omega = 2 / (2 - _iter_acceleration*_iter_acceleration);
            else
                omega = 1;

            // iterate through particles and solve system
            /** TODO: When graph coloring is used, we will iterate through all the energies in the color. */
            for (unsigned p_idx : _ctx->particles)
            {
                _solveParticle(p_idx, dt);
            }

            // Chebyshev acceleration
            for (unsigned p_idx : _ctx->particles)
            {
                _particleChebyshevAcceleration(p_idx, omega);
            }

            // for all hard constraint energies with Lagrange multipliers, update the multipliers and stiffnesses after the iteration
            _ctx->energies.forEachHardConstraintEnergyType([&] (const auto& pool) {
                for (unsigned e_idx : pool)
                {
                    /** TODO: this is hard-coded for now. Make type detectino automatic */
                    Energy::HardConstraintEnergySolver<Energy::GroundCollisionEnergyPool, Energy::GroundCollisionConstraintSolver>::updateAfterIteration(
                        e_idx,
                        pool,
                        _ctx->particles
                    );
                }
            });
        }

        // update particle velocities after iteration
        for (unsigned p_idx : _ctx->particles)
        {
            _particleVelocityUpdate(p_idx, dt);
        }
    }

private:
    void _particleInertialUpdate(unsigned p_idx, const Vec3r& a_ext, Real dt)
    {
        Vec3r p = _ctx->particles.positions[p_idx];
        Vec3r v = _ctx->particles.velocities[p_idx];
        Vec3r v_prev = _ctx->particles.previous_velocities[p_idx];

        // use adaptive initialization (Sec 3.7 in VBD paper)
        Vec3r a = (v - v_prev) / dt;

        Real a_ext_norm = a_ext.norm();
        Real a_along_a_ext = a.dot(a_ext) / a_ext_norm;

        Real a_tilde;
        if (a_along_a_ext > a_ext_norm)
            a_tilde = 1;
        else if (a_along_a_ext < 0)
            a_tilde = 0;
        else
            a_tilde = a_along_a_ext / a_ext_norm;

        Vec3r a_tilde_vec = a_tilde * a_ext;

        // compute the inertially predicted position
        _ctx->particles.inertial_positions[p_idx] = p + dt*v + dt*dt*a_ext;

        // move particle to its initialized position
        _ctx->particles.positions[p_idx] += dt*v + dt*dt*a_tilde_vec;

        // mark particle as not in collision at the beginning of the time step
        _ctx->particles.in_collision[p_idx] = false;

        // initialize the previous iteration positions
        _ctx->particles.last_iter_positions[p_idx] = _ctx->particles.positions[p_idx];
        _ctx->particles.last_last_iter_positions[p_idx] = _ctx->particles.positions[p_idx];
    }

    void _solveParticle(unsigned p_idx, Real dt)
    {
        unsigned adj_start = _ctx->adjacency.adjOffsets[p_idx];
        unsigned adj_end   = _ctx->adjacency.adjOffsets[p_idx + 1];

        Vec3r grad = Vec3r::Zero();
        Mat3r hess = Mat3r::Zero();

        // accumulate Hessians and gradients from energies
        for (unsigned e = adj_start; e < adj_end; e++) 
        {
            const ParticleAdjacency::Entry& entry = _ctx->adjacency.adj_entries[e];
            if (entry.energy_type == EnergyType::NEO_HOOKEAN)
            {
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
            else if (entry.energy_type == EnergyType::GROUND_COLLISION) 
            {
                Energy::GroundCollisionEnergySolver::accumulate(
                    entry.energy_idx,
                    _ctx->energies.ground_collision,
                    _ctx->particles,
                    entry.local_vertex_idx,
                    hess,
                    grad,
                    dt
                );
            }
        }

        // assemble LHS and RHS of single-particle system
        Real mass = _ctx->particles.masses[p_idx];
        const Vec3r& p = _ctx->particles.positions[p_idx];
        const Vec3r& y = _ctx->particles.inertial_positions[p_idx];

        Vec3r RHS = -mass / (dt*dt) * (p - y) - grad;
        Mat3r LHS = mass / (dt*dt) * Mat3r::Identity() + hess;

        Vec3r dx = LHS.partialPivLu().solve(RHS);

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