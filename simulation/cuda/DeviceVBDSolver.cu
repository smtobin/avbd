#include "simulation/cuda/DeviceVBDSolver.cuh"

namespace Sim
{

namespace
{

// Kernel definitions private to this implementation file.

__global__ void particleInertialUpdateKernel(DeviceParticlePoolView v, unsigned count, Vec3r a_ext, Real dt)
{
    unsigned p_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (p_idx >= count)
        return;

    // skip fixed particles
    if (v.fixed[p_idx])
        return;
    
    // === Positional update ===
    // assumed gravity
    Vec3r a_ext(0, -9.81, 0);

    Vec3r& p = v.positions[p_idx];
    Vec3r& y = v.inertial_positions[p_idx];
    const Vec3r& v = v.velocities[p_idx];
    const Vec3r& v_prev = v.previous_velocities[p_idx];

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

    // initialize the previous iteration positions
    v.last_iter_positions[p_idx] = p;
    v.last_last_iter_positions[p_idx] = p;

    // === Rotational update ===
    // skip if particle is not rotational
    if (!v.isOriented(p_idx))
        return;
    
    Quaternion& q = v.rotation(p_idx);
    Quaternion& q_inertial = v.inertialRotation(p_idx);
    const Vec3r& w = v.angularVelocity(p_idx);
    const Vec3r& I = v.rotationalInertia(p_idx);
    Vec3r I_inv(1.0/I[0], 1.0/I[1], 1.0/I[2]);


    /** no need to do adaptive acceleration since assuming no external angular acceleration */
    // const Vec3r& w = v.angularVelocity(p_idx);
    // const Vec3r& w_prev = v.previousAngularVelocity(p_idx);

    // Vec3r a = (w - w_prev) / dt;

    // compute inertially predicted rotation
    Vec3r alpha_ext = Vec3r(0,0,0); //Vec3r::Zero();
    Vec3r dq_inertial = dt*w + dt*dt*(alpha_ext - I_inv.cwiseProduct( ( w.cross( I.cwiseProduct(w) ) ) ) ); 
    q_inertial = Math::Plus_S3(q, dq_inertial);

    // move particle to its initialized position
    Vec3r dq = dt*w;
    q = Math::Plus_S3(q, dq);

    // initialize previous iteration rotations
    v.lastIterRotation(p_idx) = q;
    v.lastLastIterRotation(p_idx) = q;
}

__global__ void particleChebyshevAccelerationKernel(DeviceParticlePoolView v, unsigned count, Real omega)
{
    unsigned p_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (p_idx >= count)
        return;

    v.positions[p_idx] =
        omega * (v.positions[p_idx] - v.last_last_iter_positions[p_idx]) + v.last_last_iter_positions[p_idx];

    if (v.isOriented(p_idx))
    {
        Quaternion& q = v.rotation(p_idx);
        const Quaternion& last_last_q = v.lastLastIterRotation(p_idx);

        q = Math::Plus_S3(last_last_q, omega * Math::Minus_S3(q, last_last_q));
    }
    

    // update the previous iteration positions
    v.last_last_iter_positions[p_idx] = v.last_iter_positions[p_idx];
    v.last_iter_positions[p_idx] = v.positions[p_idx];

    if (v.isOriented(p_idx))
    {
        v.lastLastIterRotation(p_idx) = v.lastIterRotation(p_idx);
        v.lastIterRotation(p_idx) = v.rotation(p_idx);
    }
}

__global__ void particleVelocityUpdateKernel(DeviceParticlePoolView v, unsigned count, Real dt)
{
    unsigned p_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (p_idx >= count)
        return;

    v.previous_velocities[p_idx] = v.velocities[p_idx];
    v.velocities[p_idx] = 
        (v.positions[p_idx] - v.previous_positions[p_idx]) / dt;

    v.previous_positions[p_idx] = v.positions[p_idx];

    if (v.isOriented(p_idx))
    {
        Vec3r& prev_w = v.previousAngularVelocity(p_idx);
        Vec3r& w = v.angularVelocity(p_idx);

        const Quaternion& q = v.rotation(p_idx);
        Quaternion& prev_q = v.previousRotation(p_idx);

        prev_w = w;
        w = Math::Minus_S3(q, prev_q) / dt;

        prev_q = q;
    }
}

/** WARNING: Copied from Claude */
__global__ void solveColorKernel(SolverCtx ctx, const unsigned* __restrict__ color_particles,
                                 unsigned n, Real dt)
{
    const unsigned warp = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const unsigned lane = threadIdx.x & 31;
    if (warp >= n) return;                    // warp-uniform, safe
    const unsigned p = color_particles[warp];
    if (ctx.particles.fixed[p]) return;       // also warp-uniform, safe

    if (ctx.particles.isOriented(p))
    {
        Vec6r grad = Vec6r::Zero();
        Mat6r hess = Mat6r::Zero();

        accumulateStatic<6>(ctx, p, lane, hess, grad, dt);
        // accumulateDynamic<6>()
    }
    else
    {
        Vec3r grad = Vec3r::Zero();
        Mat3r hess = Mat3r::Zero();
        accumulateStatic<3>(ctx, p, lane, hess, grad, dt);
        // accumulateDynamic<3>()

        warpReduceSym<3>(hess, grad);           // every lane now holds the totals
        if (lane == 0)
            assembleSolveWrite<3>(ctx, p, hess, grad, dt);   // your existing tail
    }

    /** TODO: (10/26/26) write/understand this function */
    

    
}

/** WARNING: Copied from Claude */
template<int DOF>
__device__ void accumulateStatic(const SolverCtx& ctx, unsigned p, unsigned lane,
                                 Mat<DOF>& hess, Vec<DOF>& grad, Real dt)
{
    const auto& offs = ctx.static_adjacency.e_offsets[p];
    const unsigned adj_end = ctx.static_adjacency.e_offsets[p+1][0];

    Energy::ForEachStaticEnergy([&]<StaticEnergyType E>() {
        constexpr unsigned t = static_cast<unsigned>(E);
        using Solver = SolverFor<E>::type;
        if constexpr ((DOF == 6 && Solver::SupportsOriented) || (DOF == 3 && Solver::SupportsPositional))
        {
            const unsigned begin = offs[t];
            const unsigned end   = (t + 1 == kNumStatic) ? adj_end : offs[t + 1];
            const auto& pool = ctx.energies.template get<E>();

            for (unsigned e = begin + lane; e < end; e += 32)
            {
                const auto& entry = ctx.static_adjacency.e_entries[e];
                Solver::template accumulate<DOF>(entry.energy_idx, pool, ctx.particles,
                                                 entry.local_vertex_idx, hess, grad, dt);
            }
        }
    });
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