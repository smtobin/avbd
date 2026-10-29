#include "simulation/cuda/DeviceVBDSolver.cuh"
#include "common/EnergyUtils.hpp"

#include "energy/elastic/NeoHookeanEnergySolver.hpp"
#include "energy/elastic/CosseratRodEnergySolver.hpp"
#include "energy/collision/GroundCollisionEnergySolver.hpp"
#include "energy/collision/RigidBodyGroundCollisionEnergySolver.hpp"
#include "energy/collision/TriangleRigidCollisionEnergySolver.hpp"
#include "energy/collision/TriangleRodCollisionEnergySolver.hpp"
#include "energy/joint/OneSidedFixedJointEnergySolver.hpp"

namespace Sim
{

namespace
{

// Kernel definitions private to this implementation file.

template <typename PoolView>
__global__ void updateAfterTimeStepKernel(PoolView view, ParticlePoolView particles)
{
    unsigned e_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (e_idx >= view.counters->highest_index+1)
        return;

    if (!view.active[e_idx])
        return;

    if constexpr (Energy::HasUpdateContactPoints<typename PoolView::PoolType::SolverType>)
    {
        PoolView::PoolType::SolverType::updateContactPoints(e_idx, view, particles);
    }
    PoolView::PoolType::SolverType::updateAfterTimeStep(e_idx, view, particles);
}

template <typename PoolView>
__global__ void updateAfterIterationKernel(PoolView view, ParticlePoolView particles)
{
    unsigned e_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (e_idx >= view.counters->highest_index+1)
        return;

    if (!view.active[e_idx])
        return;

    PoolView::PoolType::SolverType::updateAfterIteration(e_idx, view, particles);
}

__global__ void particleInertialUpdateKernel(ParticlePoolView view, unsigned count, Vec3r a_ext, Real dt)
{
    unsigned p_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (p_idx >= count)
        return;

    // skip fixed particles
    if (view.fixed[p_idx])
        return;
    
    // === Positional update ===
    Vec3r& p = view.positions[p_idx];
    Vec3r& y = view.inertial_positions[p_idx];
    const Vec3r& v = view.velocities[p_idx];
    const Vec3r& v_prev = view.previous_velocities[p_idx];

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
    view.last_iter_positions[p_idx] = p;
    view.last_last_iter_positions[p_idx] = p;

    // === Rotational update ===
    // skip if particle is not rotational
    if (!view.isOriented(p_idx))
        return;
    
    Quaternion& q = view.rotation(p_idx);
    Quaternion& q_inertial = view.inertialRotation(p_idx);
    const Vec3r& w = view.angularVelocity(p_idx);
    const Vec3r& I = view.rotationalInertia(p_idx);
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
    view.lastIterRotation(p_idx) = q;
    view.lastLastIterRotation(p_idx) = q;
}

__global__ void particleChebyshevAccelerationKernel(ParticlePoolView v, unsigned count, Real omega)
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

__global__ void particleVelocityUpdateKernel(ParticlePoolView v, unsigned count, Real dt)
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


template<int DOF>
__device__ void accumulateStatic(
    const SimulationContextView& ctx,
    unsigned p,
    unsigned lane,
    MatNr<DOF>& hess,
    VecNr<DOF>& grad,
    Real dt
)
{
    const auto& offs = ctx.static_adjacency.e_offsets[p];
    const unsigned adj_end = ctx.static_adjacency.e_offsets[p+1][0];
    unsigned num_energies = static_cast<unsigned>(StaticEnergyType::count);

    Energy::ForEachStaticEnergy([&]<StaticEnergyType E>() {
        constexpr unsigned t = static_cast<unsigned>(E);
        using Solver = SolverFor<E>::type;
        if constexpr ((DOF == 6 && Solver::SupportsOriented) || (DOF == 3 && Solver::SupportsPositional))
        {
            const unsigned begin = offs[t];
            const unsigned end   = (t + 1 == num_energies) ? adj_end : offs[t + 1];
            const auto& pool_view = ctx.energies.template get<E>();

            for (unsigned e = begin + lane; e < end; e += 32)
            {
                const auto& entry = ctx.static_adjacency.e_entries[e];
                Solver::template accumulate<DOF>(entry.energy_idx, pool_view, ctx.particles,
                                                 entry.local_vertex_idx, hess, grad, dt);
            }
        }
    });
}

__device__ void assembleSolveWrite(
    SimulationContextView& ctx,
    unsigned p_idx, 
    const Mat3r& hess, 
    const Vec3r& grad, 
    Real dt
)
{
    Real mass = ctx.particles.masses[p_idx];
    const Vec3r& p = ctx.particles.positions[p_idx];
    const Vec3r& y = ctx.particles.inertial_positions[p_idx];

    Vec3r RHS = -mass / (dt*dt) * (p - y) - grad;
    Mat3r LHS = mass / (dt*dt) * Mat3r::Identity() + hess;

    

    Vec3r dx = Vec3r::Zero();
    bool success = LHS.solve(RHS, dx);
    /** TODO: (10/29/26) Set error flags when solver fails. */
    // if (!success)
    //     throw std::runtime_error("Failed to solve 3x3 particle system!");
        
    // std::cout << "Particle " << p_idx << ": " << "\n  LHS:\n" << LHS << "\n  RHS: " << RHS.transpose() << "\n  dx: " << dx.transpose() << std::endl;
    // Vec3r dx = LHS.inverse() * RHS;

    // std::cout << "dx: " << dx.transpose() << std::endl;

    /** TODO: (10/29/26) Handle intra-color conflicts. */
    // if this particle has an intra-color conflict, buffer its update
    // if (ctx.coloring.is_conflicted[p_idx])
    //     ctx.particles.buffered_positions[p_idx] = _ctx.particles.positions[p_idx] + dx;
    // otherwise handle it normally
    // else
        ctx.particles.positions[p_idx] += dx;
}

__device__ void assembleSolveWriteOriented(
    SimulationContextView& ctx,
    unsigned p_idx,
    const Mat6r& hess,
    const Vec6r& grad,
    Real dt
)
{
    // extract quantities from pool
    Real mass = ctx.particles.masses[p_idx];
    const Vec3r& I = ctx.particles.rotationalInertia(p_idx);
    Vec3r& p = ctx.particles.positions[p_idx];
    const Vec3r& p_inertial = ctx.particles.inertial_positions[p_idx];
    Quaternion& q = ctx.particles.rotation(p_idx);
    const Quaternion& q_inertial = ctx.particles.inertialRotation(p_idx);

    // form RHS
    Vec6r RHS;
    RHS.head<3>() = -1/(dt*dt) * mass * (p - p_inertial);
    RHS.tail<3>() = -1/(dt*dt) * I.cwiseProduct( Math::Minus_S3(q, q_inertial) );
    RHS -= grad;

    // form LHS
    Mat6r LHS = Mat6r::Zero();
    Vec3r diag1( mass / (dt*dt) );
    LHS(0,0) = diag1[0]; LHS(1,1) = diag1[1]; LHS(2,2) = diag1[2];
    Vec3r diag2 = I / (dt*dt);
    LHS(3,3) = diag2[0]; LHS(4,4) = diag2[1]; LHS(5,5) = diag2[2];
    LHS += hess;

    Vec6r dx = Vec6r::Zero();
    bool success = LHS.solve(RHS, dx);
    /** TODO: (10/29/26) Set error flags if solver fails. */
    // if (!success)
    //     throw std::runtime_error("Failed to solve 6x6 particle system!");

    Vec3r dp = dx.head<3>();
    Vec3r dq = dx.tail<3>();

    /** TODO: (10/29/26) Handle intra-color conflicts. */
    // if this particle has an intra-color conflict, buffer its update
    // if (ctx.coloring.is_conflicted[p_idx])
    // {
    //     ctx.particles.buffered_positions[p_idx] = p + dp;
    //     ctx.particles.bufferedRotation(p_idx) = Math::Plus_S3(q, dq);
    // }
    // // otherwise update normally
    // else
    // {
        p += dp;
        q = Math::Plus_S3(q, dq);
    // }
}

template<typename T>
__device__ __forceinline__ T warpSum(T v) {
    #pragma unroll
    for (int o = 16; o > 0; o >>= 1)
        v += __shfl_xor_sync(0xffffffffu, v, o);
    return v;
}

template<int DOF>
__device__ __forceinline__ void warpReduceSym(MatNr<DOF>& H, VecNr<DOF>& g)
{
    #pragma unroll
    for (int i = 0; i < DOF; ++i) {
        #pragma unroll
        for (int j = i; j < DOF; ++j) H(i, j) = warpSum(H(i, j));
        g[i] = warpSum(g[i]);
    }
    #pragma unroll
    for (int i = 0; i < DOF; ++i)
        #pragma unroll
        for (int j = 0; j < i; ++j) H(i, j) = H(j, i);
}


__global__ void solveColorKernel(
    SimulationContextView ctx,
    unsigned c_idx,
    Real dt
)
{
    const unsigned color_start = ctx.coloring.color_offsets[c_idx];

    const unsigned warp = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const unsigned lane = threadIdx.x & 31;
    if (warp >= ctx.coloring.color_counts[c_idx]) return;                    // warp-uniform, safe
    const unsigned p = ctx.coloring.work_list[color_start + warp];
    if (ctx.particles.fixed[p]) return;       // also warp-uniform, safe

    if (ctx.particles.isOriented(p))
    {
        Vec6r grad = Vec6r::Zero();
        Mat6r hess = Mat6r::Zero();
        accumulateStatic<6>(ctx, p, lane, hess, grad, dt);
        // accumulateDynamic<6>()

        warpReduceSym<6>(hess, grad);
        if (lane == 0)
            assembleSolveWriteOriented(ctx, p, hess, grad, dt);
    }
    else
    {
        Vec3r grad = Vec3r::Zero();
        Mat3r hess = Mat3r::Zero();
        accumulateStatic<3>(ctx, p, lane, hess, grad, dt);
        // accumulateDynamic<3>()

        warpReduceSym<3>(hess, grad);           // every lane now holds the totals
        if (lane == 0)
            assembleSolveWrite(ctx, p, hess, grad, dt);   // your existing tail
    }
    
}
} // anonymous namespace

DeviceVBDSolver::DeviceVBDSolver(
    DeviceSimulationContext* ctx,
    unsigned solver_iters,
    Real iter_acceleration)
    : _ctx(ctx->view())
    , _params(&ctx->params)
    , _solver_iters(solver_iters)
    , _iter_acceleration(iter_acceleration)
{
}

void DeviceVBDSolver::solve()
{
    int num_particles = _ctx.particles.counters->highest_index + 1;

    // update all energies after time step
    _ctx.energies.forEachEnergyType([&] (auto& view) {
        int threads = 256;
        /** TODO: (10/28/26) for now, using the pool capacity as upper bound on energies */
        int e_blocks = (view.capacity + threads - 1) / threads;
        updateAfterTimeStepKernel<<<e_blocks, threads>>>(view, _ctx.particles);
    });
    
    
    // compute inertial update and initialization for each particle
    int threads = 256;
    int p_blocks = (num_particles + threads - 1) / threads;
    particleInertialUpdateKernel<<<p_blocks, threads>>>(_ctx.particles, num_particles, Vec3r(0,-9.81,0), _params->dt);

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
        for (unsigned c = 0; c < _ctx.num_colors; c++)
        {
            /** TODO: (10/28/26) for now, using upper bound on particles in this color */
            solveColorKernel<<<p_blocks, threads>>>(_ctx, c, _params->dt);
        }

        // Chebyshev acceleration
        particleChebyshevAccelerationKernel<<<p_blocks, threads>>>(_ctx.particles, num_particles, omega);

        // update Lagrange multipliers
        _ctx.energies.forEachHardConstraintEnergyType([&] (auto& view) {
            int threads = 256;
            /** TODO: (10/28/26) for now, using the pool capacity as upper bound on energies */
            int e_blocks = (view.capacity + threads - 1) / threads;
            updateAfterIterationKernel<<<e_blocks, threads>>>(view, _ctx.particles);
        });
    }

    // velocity update
    particleVelocityUpdateKernel<<<p_blocks, threads>>>(_ctx.particles, num_particles, _params->dt);
}

} // namespace Sim