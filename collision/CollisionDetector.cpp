#include "collision/CollisionDetector.hpp"

#include "simulation/SimulationContext.hpp"
#include "common/WorkerThreadContext.hpp"

#include "energy/collision/TriangleRigidCollisionEnergySolver.hpp"
#include "energy/collision/TriangleRodCollisionEnergySolver.hpp"

#include "common/Math.hpp"
#include "common/Algorithm.hpp"

namespace Collision
{

CollisionDetector::CollisionDetector(unsigned capacity)
{
    _potential_collisions.reserve(capacity);
    _prev_detected_collisions.reserve(capacity);
    _prev_sorted_order.reserve(capacity);
    _cur_detected_collisions.reserve(capacity);
    _cur_sorted_order.reserve(capacity);
}

void CollisionDetector::detectCollisionsAndRecolor(Sim::SimulationContext& ctx)
{
    // move the detected collision from last frame to the previous frame
    _prev_detected_collisions.swap(_cur_detected_collisions);
    _prev_sorted_order.swap(_cur_sorted_order);
    // and clear the current detected collisions
    _cur_detected_collisions.clear();
    _cur_sorted_order.clear();
    // clear the potential collisions
    _potential_collisions.clear();


    // build BVH
    _lbvh_builder.buildBVH(ctx.particles, ctx.collision_pool, ctx.lbvh, ctx.params.dt);

    // traverse BVH for potential collisions
    _lbvh_traversal.traverseSelfIterative(ctx.lbvh, ctx.lbvh.root, ctx.lbvh.root, _potential_collisions);

    // narrow-phase collision detection
    _narrowPhaseCollisionDetection(ctx);

    // handle detected collisions
    // sort the current detected collisions by key
    // merge lists, and create/destroy collision constraints accordingly
    _handleDetectedCollisions(ctx);

    // rebuild adjacency
    /** TODO: (07/22/26) Do this incrementally? Handle collision constraints separately? Anything to not recompute adjacency from scratch each time. */
    // ctx.static_adjacency.buildAdjacency(ctx.particles, ctx.energies);
    ctx.dynamic_adjacency.buildAdjacency(ctx.particles, ctx.energies);
    // ctx.coloring.buildColorList(ctx.adjacency, ctx.particles.totalSize());
}

void CollisionDetector::detectCollisionsAndRecolor_Parallel(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, Sim::SimulationContext& ctx)
{
    if (w_ctx.idx == 0)
    {
        // move the detected collision from last frame to the previous frame
        _prev_detected_collisions.swap(_cur_detected_collisions);
        _prev_sorted_order.swap(_cur_sorted_order);
        // and clear the current detected collisions
        _cur_detected_collisions.clear();
        _cur_sorted_order.clear();
        // clear the potential collisions
        _potential_collisions.clear();
    }
    w_ctx.barrier->arrive_and_wait();

    // build BVH
    _lbvh_builder.buildBVH_Parallel(w_ctx, all_worker_contexts, ctx);

    // traverse BVH for potential collisions
    _lbvh_traversal.collisionBroadPhase(w_ctx, ctx.lbvh);
    _lbvh_traversal.mergeLeafPairs(w_ctx, all_worker_contexts, _potential_collisions);

    // narrow-phase collision detection
    _narrowPhaseCollisionDetection_Parallel(w_ctx, ctx);
    _mergeDetectedCollisions(w_ctx, all_worker_contexts, _cur_detected_collisions);

    if (w_ctx.idx == 0)
    {
        // handle detected collisions
        // sort the current detected collisions by key
        // merge lists, and create/destroy collision constraints accordingly
        _handleDetectedCollisions(ctx);

        
    }
    w_ctx.barrier->arrive_and_wait();

    // ctx.adjacency.buildAdjacency_Parallel(w_ctx, all_worker_contexts, ctx.particles, ctx.energies);
    if (w_ctx.idx == 0)
    {
        // rebuild adjacency
        // std::cout << "=== Adjacency parallel ===" << std::endl;
        // for (unsigned v : ctx.particles)
        // {
        //     std::cout << "Particle " << v << ": " << std::endl;
        //     unsigned start = ctx.adjacency.e_offsets[v][0];
        //     unsigned end = ctx.adjacency.e_offsets[v+1][0];
        //     for (unsigned i = start; i < end; i++)
        //     {
        //         std::cout << " Entry " << i << ": " << "e_idx=" << ctx.adjacency.e_entries[i].energy_idx << " type=" << static_cast<unsigned>(ctx.adjacency.e_entries[i].energy_type) << std::endl;
        //     }
        // }
        /** TODO: (07/22/26) Do this incrementally? Handle collision constraints separately? Anything to not recompute adjacency from scratch each time. */
        // ctx.static_adjacency.buildAdjacency(ctx.particles, ctx.energies);
        ctx.dynamic_adjacency.buildAdjacency(ctx.particles, ctx.energies);
        ctx.coloring.incrementalRecoloring(ctx.static_adjacency, ctx.dynamic_adjacency, ctx.particles.totalSize());

        // std::cout << "=== Adjacency serial ===" << std::endl;
        // for (unsigned v : ctx.particles)
        // {
        //     std::cout << "Particle " << v << ": " << std::endl;
        //     unsigned start = ctx.adjacency.e_offsets[v][0];
        //     unsigned end = ctx.adjacency.e_offsets[v+1][0];
        //     for (unsigned i = start; i < end; i++)
        //     {
        //         std::cout << " Entry " << i << ": " << "e_idx=" << ctx.adjacency.e_entries[i].energy_idx << " type=" << static_cast<unsigned>(ctx.adjacency.e_entries[i].energy_type) << std::endl;
        //     }
        // }

        // ctx.coloring.buildInitialColorList(ctx.static_adjacency, ctx.dynamic_adjacency, ctx.particles.totalSize());
    }
    w_ctx.barrier->arrive_and_wait();
}


bool CollisionDetector::_shouldSkip(const CollisionPrimitivePool& cpool, unsigned pi, unsigned pj)
{
    // always test different objects
    if (cpool.object_id[pi] != cpool.object_id[pj])
        return false;

    // check if primitives share vertex
    // if so, skip the collision
    const auto& inds_i = cpool.particle_indices[pi];
    const auto& inds_j = cpool.particle_indices[pj];
    unsigned cnt_i = cpool.num_particles[pi];
    unsigned cnt_j = cpool.num_particles[pj];

    // std::cout << " Checking primitives " << pi << " and " << pj << "..." << std::endl;
    // std::cout << "  Primitive i: ";
    // for (unsigned ind : inds_i) std::cout << ind << ", ";
    // std::cout << std::endl << "  Primitive j: ";
    // for (unsigned ind : inds_j) std::cout << ind << ", ";
    // std::cout << std::endl;
    // std::cout << "  AABB i: " << cpool.aabb[pi] << std::endl;
    // std::cout << "  AABB j: " << cpool.aabb[pj] << std::endl;

    for (unsigned i = 0; i < cnt_i; i++)
    {
        for (unsigned j = 0; j < cnt_j; j++)
        {
            if (inds_i[i] == inds_j[j])
                return true;
        }
    }

    return false;
}

void CollisionDetector::_processPotentialCollision(Sim::SimulationContext& ctx, unsigned a, unsigned b, std::vector<DetectedCollision>& detected_collisions)
{
    if (_shouldSkip(ctx.collision_pool, a, b))
            return;

        CollisionGeometryType type_a = ctx.collision_pool.getCollisionGeometryType(a);
        CollisionGeometryType type_b = ctx.collision_pool.getCollisionGeometryType(b);

        if (type_a > type_b)
        {
            std::swap(type_a, type_b);
            std::swap(a, b);
        }

        switch (_makeCollisionKey(type_a, type_b))
        {
            
            case _makeCollisionKey(CollisionGeometryType::Triangle, CollisionGeometryType::Triangle):
            {
                _triangleTriangle(ctx, a, b, detected_collisions);
                break;
            }
            case _makeCollisionKey(CollisionGeometryType::Triangle, CollisionGeometryType::Sphere):
            {
                _triangleSphere(ctx, a, b, detected_collisions);
                break;
            }
            case _makeCollisionKey(CollisionGeometryType::Triangle, CollisionGeometryType::RodSegment):
            {
                _triangleRod(ctx, a, b, detected_collisions);
                break;
            }
            case _makeCollisionKey(CollisionGeometryType::RodSegment, CollisionGeometryType::RodSegment):
            {
                _rodRod(ctx, a, b, detected_collisions);
                break;
            }
            case _makeCollisionKey(CollisionGeometryType::Sphere, CollisionGeometryType::Sphere):
            {
                _sphereSphere(ctx, a, b, detected_collisions);
                break;
            }

            default:
            {
                throw std::runtime_error("Unsupported collision!");
            }
            
        }
}

void CollisionDetector::_narrowPhaseCollisionDetection(Sim::SimulationContext& ctx)
{
    for (const auto& potential_collision : _potential_collisions)
    {
        // potential collisions are pairs of LBVH leaf node indices, which correspond to the collision pool's sorted order
        // extract the primitive indices by indexing in the sorted order array
        unsigned a = ctx.collision_pool.sorted_order[potential_collision.first];
        unsigned b = ctx.collision_pool.sorted_order[potential_collision.second];

        _processPotentialCollision(ctx, a, b, _cur_detected_collisions);
    }
}

void CollisionDetector::_narrowPhaseCollisionDetection_Parallel(WorkerThreadContext& w_ctx, Sim::SimulationContext& ctx)
{
    auto [start, end] = w_ctx.computeStartEnd(_potential_collisions.size());
    w_ctx.detected_collisions.clear();

    for (unsigned c_idx = start; c_idx < end; c_idx++)
    {
        auto& potential_collision = _potential_collisions[c_idx];
        // potential collisions are pairs of LBVH leaf node indices, which correspond to the collision pool's sorted order
        // extract the primitive indices by indexing in the sorted order array
        unsigned a = ctx.collision_pool.sorted_order[potential_collision.first];
        unsigned b = ctx.collision_pool.sorted_order[potential_collision.second];

        _processPotentialCollision(ctx, a, b, w_ctx.detected_collisions);
    }
    w_ctx.barrier->arrive_and_wait();
}

void CollisionDetector::_mergeDetectedCollisions(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, std::vector<DetectedCollision>& merged_detected_collisions)
{
    // main thread compute per-thread offsets
    if (w_ctx.idx == 0)
    {
        _merge_offsets.resize(WorkerThreadContext::NUM_THREADS + 1);
        _merge_offsets[0] = 0;
        for (unsigned t = 0; t < WorkerThreadContext::NUM_THREADS; t++)
        {
            _merge_offsets[t+1] = _merge_offsets[t] + all_worker_contexts[t].detected_collisions.size();
        }
        merged_detected_collisions.resize(_merge_offsets.back());
    }
    w_ctx.barrier->arrive_and_wait();

    // each thread copies its own buffer into its assigned slice
    std::copy(
        w_ctx.detected_collisions.begin(), w_ctx.detected_collisions.end(),
        merged_detected_collisions.begin() + _merge_offsets[w_ctx.idx]
    );
    w_ctx.barrier->arrive_and_wait();
}


void CollisionDetector::_handleDetectedCollisions(Sim::SimulationContext& ctx)
{
    // sort detected collision list by key
    Algorithm::radixSort(_cur_detected_collisions, _cur_sorted_order, _cur_detected_collisions.size(), 
        [](const DetectedCollision& c) { return c.key; } );

    const unsigned cur_size  = static_cast<unsigned>(_cur_detected_collisions.size());
    const unsigned prev_size = static_cast<unsigned>(_prev_detected_collisions.size());

    unsigned cur_idx = 0;
    unsigned prev_idx = 0;

    while (cur_idx < cur_size && prev_idx < prev_size)
    {
        DetectedCollision& cur  = _cur_detected_collisions[_cur_sorted_order[cur_idx]];
        DetectedCollision& prev = _prev_detected_collisions[_prev_sorted_order[prev_idx]];

        if (cur.key < prev.key)
        {
            // no previous collision with this key - add a new collision
            _addCollision(ctx, cur);
            cur_idx++;
        }
        else if (cur.key > prev.key)
        {
            // previous collision no longer detected this frame - remove it
            _removeCollision(ctx, prev);
            prev_idx++;
        }
        else // cur.key == prev.key
        {
            if (cur.gen1 == prev.gen1 && cur.gen2 == prev.gen2)
            {
                // same collision was in previous frame - update normal, contact points, etc.
                cur.e_idx = prev.e_idx;
                _updateCollision(ctx, cur);
            }
            else
            {
                // same key, but a stale slot got reused - treat as unrelated
                _removeCollision(ctx, prev);
                _addCollision(ctx, cur);
            }
            cur_idx++;
            prev_idx++;
        }
    }

    // anything left in cur has no counterpart - all new
    for (; cur_idx < cur_size; cur_idx++)
        _addCollision(ctx, _cur_detected_collisions[_cur_sorted_order[cur_idx]]);

    // anything left in prev has no counterpart - all removed
    for (; prev_idx < prev_size; prev_idx++)
        _removeCollision(ctx, _prev_detected_collisions[_prev_sorted_order[prev_idx]]);
}

void CollisionDetector::_addCollision(Sim::SimulationContext& ctx, DetectedCollision& collision)
{
    switch (collision.type)
    {
        case DetectedCollisionType::TriangleTriangle_VertexFace:
        {
            throw std::runtime_error("TriangleTriangle Vertex-Face collisions not implemented!");
            break;
        }
        case DetectedCollisionType::TriangleTriangle_EdgeEdge:
        {
            throw std::runtime_error("TriangleTriangle Edge-Edge collisions not implemented!");
            break;
        }
        case DetectedCollisionType::TriangleRigid:
        {
            /** TODO: (08/05/26) set coefficients of friction based on material properties */

            // set starting collision stiffness based on minimum particle mass involved
            Real min_mass = std::min({
                ctx.particles.masses[collision.TriangleRigid.tri[0]],
                ctx.particles.masses[collision.TriangleRigid.tri[1]],
                ctx.particles.masses[collision.TriangleRigid.tri[2]],
                ctx.particles.masses[collision.TriangleRigid.rb]
            });
            Real k_start = min_mass / (ctx.params.dt * ctx.params.dt);

            Vec3r t, b;
            Math::completeOrthonormalBasisGivenNormal(collision.normal, t, b);
            unsigned slot = ctx.energies.triangle_rigid_collision.addEnergy(
                collision.TriangleRigid.tri[0],
                collision.TriangleRigid.tri[1],
                collision.TriangleRigid.tri[2],
                collision.TriangleRigid.rb,
                nullptr,
                collision.normal,
                t,
                b,
                collision.TriangleRigid.barys,
                collision.TriangleRigid.cp_rb_local,
                k_start,
                0.4, 
                0.2
            );
            collision.e_idx = slot;

            // add particles involved in this collision to the coloring dirty list for recoloring
            ctx.coloring.markParticleDirty(collision.TriangleRigid.tri[0]);
            ctx.coloring.markParticleDirty(collision.TriangleRigid.tri[1]);
            ctx.coloring.markParticleDirty(collision.TriangleRigid.tri[2]);
            ctx.coloring.markParticleDirty(collision.TriangleRigid.rb);
            break;
        }
        case DetectedCollisionType::TriangleRod:
        {
            /** TODO: (08/26/26) set coefficients of friction based on material properties */

            // set starting collision stiffness based on minimum particle mass involved
            Real min_mass = std::min({
                ctx.particles.masses[collision.TriangleRod.tri[0]],
                ctx.particles.masses[collision.TriangleRod.tri[1]],
                ctx.particles.masses[collision.TriangleRod.tri[2]],
                ctx.particles.masses[collision.TriangleRod.rod[0]],
                ctx.particles.masses[collision.TriangleRod.rod[1]]
            });
            Real k_start = min_mass / (ctx.params.dt * ctx.params.dt);

            Vec3r t, b;
            Math::completeOrthonormalBasisGivenNormal(collision.normal, t, b);
            unsigned slot = ctx.energies.triangle_rod_collision.addEnergy(
                collision.TriangleRod.tri[0],
                collision.TriangleRod.tri[1],
                collision.TriangleRod.tri[2],
                collision.TriangleRod.rod[0],
                collision.TriangleRod.rod[1],
                collision.normal,
                t,
                b,
                collision.TriangleRod.barys,
                collision.TriangleRod.s,
                collision.TriangleRod.cp_rod_local,
                k_start,
                0.0, 
                0.0
            );
            collision.e_idx = slot;

            // add particles involved in this collision to the coloring dirty list for recoloring
            ctx.coloring.markParticleDirty(collision.TriangleRod.tri[0]);
            ctx.coloring.markParticleDirty(collision.TriangleRod.tri[1]);
            ctx.coloring.markParticleDirty(collision.TriangleRod.tri[2]);
            ctx.coloring.markParticleDirty(collision.TriangleRod.rod[0]);
            ctx.coloring.markParticleDirty(collision.TriangleRod.rod[1]);
            break;
        }
        default:
        {
            throw std::runtime_error("CollisionDetector::_addCollision - unsupported DetectedCollisionType!");
            break;
        }
    }
}

void CollisionDetector::_removeCollision(Sim::SimulationContext& ctx, DetectedCollision& collision)
{
    switch (collision.type)
    {
        case DetectedCollisionType::TriangleTriangle_VertexFace:
        {
            throw std::runtime_error("TriangleTriangle Vertex-Face collisions not implemented!");
            break;
        }
        case DetectedCollisionType::TriangleTriangle_EdgeEdge:
        {
            throw std::runtime_error("TriangleTriangle Edge-Edge collisions not implemented!");
            break;
        }
        case DetectedCollisionType::TriangleRigid:
        {
            ctx.energies.triangle_rigid_collision.removeEnergy(collision.e_idx);
            break;
        }
        case DetectedCollisionType::TriangleRod:
        {
            ctx.energies.triangle_rod_collision.removeEnergy(collision.e_idx);
            break;
        }
        default:
        {
            throw std::runtime_error("CollisionDetector::_removeCollision - unsupported DetectedCollisionType!");
            break;
        }
    }
}

void CollisionDetector::_updateCollision(Sim::SimulationContext& ctx, DetectedCollision& collision)
{
    switch (collision.type)
    {
        case DetectedCollisionType::TriangleTriangle_VertexFace:
        {
            throw std::runtime_error("TriangleTriangle Vertex-Face collisions not implemented!");
            break;
        }
        case DetectedCollisionType::TriangleTriangle_EdgeEdge:
        {
            throw std::runtime_error("TriangleTriangle Edge-Edge collisions not implemented!");
            break;
        }
        case DetectedCollisionType::TriangleRigid:
        {
            Energy::TriangleRigidCollisionEnergyInfo& info = ctx.energies.triangle_rigid_collision.data[collision.e_idx];
            info.cp_rb_local = collision.TriangleRigid.cp_rb_local;
            info.barys = collision.TriangleRigid.barys;
            Energy::TriangleRigidCollisionEnergySolver::updateCollisionFrame(collision.e_idx, ctx.energies.triangle_rigid_collision, collision.normal);
            break;
        }
        case DetectedCollisionType::TriangleRod:
        {
            Energy::TriangleRodCollisionEnergyInfo& info = ctx.energies.triangle_rod_collision.data[collision.e_idx];
            info.s = collision.TriangleRod.s;
            info.cp_rod_local = collision.TriangleRod.cp_rod_local;
            info.barys = collision.TriangleRod.barys;
            Energy::TriangleRodCollisionEnergySolver::updateCollisionFrame(collision.e_idx, ctx.energies.triangle_rod_collision, collision.normal);
            break;
        }
        default:
        {
            throw std::runtime_error("CollisionDetector::_updateCollision - unsupported DetectedCollisionType!");
            break;
        }
    }
}


void CollisionDetector::_triangleSphere(Sim::SimulationContext& ctx, unsigned triangle, unsigned sphere, std::vector<DetectedCollision>& detected_collisions)
{
    // std::cout << "Testing sphere-triangle collision..." << std::endl;
    const auto& triangle_idx = ctx.collision_pool.particle_indices[triangle];

    unsigned sdf_idx = ctx.collision_pool.particle_indices[sphere][0];
    unsigned sphere_idx = ctx.collision_pool.sdf_pool.particles[sdf_idx];
    const auto& sphere_params = ctx.collision_pool.sdf_pool.params[sdf_idx];

    // extract current triangle vertex positions
    const Vec3r& v1 = ctx.particles.positions[triangle_idx[0]];
    const Vec3r& v2 = ctx.particles.positions[triangle_idx[1]];
    const Vec3r& v3 = ctx.particles.positions[triangle_idx[2]];

    // extract current sphere center
    const Vec3r& p = ctx.particles.positions[sphere_idx];

    // closest point on triangle to sphere center
    Vec3r tri_cp = Math::closestPoint_PointTriangle(p, v1, v2, v3);

    Vec3r normal, cp_barys, cp_rb_local;
    if (_triangleSDF_CCD(ctx, triangle, sphere, ctx.params.dt, normal, cp_barys, cp_rb_local))
    {
        // std::cout << "Sphere-triangle collision detected!" << std::endl;
        // std::cout << "  Barys: " << cp_barys.transpose() << std::endl;
        // std::cout << "  CP on triangle: " << (v1*cp_barys[0] + v2*cp_barys[1] + v3*cp_barys[2]).transpose() << std::endl;
        // std::cout << "  CP on rigid body: " << cp_rb_local.transpose() << std::endl;
        // std::cout << "  Normal: " << normal.transpose() << std::endl;
        DetectedCollision collision{};
        collision.type = DetectedCollisionType::TriangleRigid;
        collision.key = DetectedCollision::generateKey(
            DetectedCollisionType::TriangleRigid,
            triangle,
            sphere,
            0
        );
        collision.gen1 = ctx.collision_pool.generation[triangle];
        collision.gen2 = ctx.collision_pool.generation[sphere];
        collision.normal = normal;
        collision.TriangleRigid.tri = Vec3u(triangle_idx[0], triangle_idx[1], triangle_idx[2]);
        collision.TriangleRigid.barys = cp_barys;
        collision.TriangleRigid.rb = sphere_idx;
        collision.TriangleRigid.cp_rb_local = cp_rb_local;

        detected_collisions.push_back(std::move(collision));
    }

    // for now, expand sphere radius by relative velocity
    // const Vec3r& sphere_vel = ctx.particles.velocities[sphere_idx];
    // const Vec3r& vel1 = ctx.particles.velocities[triangle_idx[0]];
    // const Vec3r& vel2 = ctx.particles.velocities[triangle_idx[1]];
    // const Vec3r& vel3 = ctx.particles.velocities[triangle_idx[2]];
    // const Vec3r avg_vel = (vel1 + vel2 + vel3) / 3;
    // const Vec3r rel_vel = sphere_vel - avg_vel;
    // Real travel = rel_vel.norm() * ctx.params.dt;

    // // check distance between closest point and sphere center
    // Vec3r diff = (tri_cp - p);
    // Real dist = diff.norm();
    // if (dist < sphere_params.sphere.radius + travel)
    // {
    //     // compute collision normal
    //     Vec3r normal;
    //     if (dist > Real(1e-6))
    //         normal = diff/dist;
    //     else
    //         normal = Vec3r(1,0,0);  // arbitrary direction

    //     // barycentric coordinates of triangle collision point
    //     Vec3r barys = Math::barycentricCoordinates(tri_cp, v1, v2, v3);

    //     // sphere conact point in local frame
    //     Vec3r sphere_cp_local = ctx.particles.rotation(sphere_idx).inverse() * normal*sphere_params.sphere.radius;

    //     DetectedCollision collision{};
    //     collision.type = DetectedCollisionType::TriangleRigid;
    //     collision.key = DetectedCollision::generateKey(
    //         DetectedCollisionType::TriangleRigid,
    //         triangle,
    //         sphere,
    //         0
    //     );
    //     collision.gen1 = ctx.collision_pool.generation[triangle];
    //     collision.gen2 = ctx.collision_pool.generation[sphere];
    //     collision.normal = normal;
    //     collision.TriangleRigid.tri = Vec3u(triangle_idx[0], triangle_idx[1], triangle_idx[2]);
    //     collision.TriangleRigid.barys = barys;
    //     collision.TriangleRigid.rb = sphere_idx;
    //     collision.TriangleRigid.cp_rb_local = sphere_cp_local;

    //     _cur_detected_collisions.push_back(std::move(collision));

    //     std::cout << "Triangle-sphere collision!" << std::endl;
    //     std::cout << "  Normal: " << normal.transpose() << std::endl;
    // }

}

void CollisionDetector::_triangleTriangle(Sim::SimulationContext& ctx, unsigned triangle1, unsigned triangle2, std::vector<DetectedCollision>& detected_collisions)
{
    /** TODO: (07/20/26) triangle-triangle collision detection */
    // std::cout << "Testing triangle-triangle collision..." << std::endl;
    // throw std::runtime_error("Triangle-triangle collision detection not implemented.");
}

void CollisionDetector::_triangleRod(Sim::SimulationContext& ctx, unsigned triangle, unsigned rod, std::vector<DetectedCollision>& detected_collisions)
{
    const auto& triangle_idx = ctx.collision_pool.particle_indices[triangle];
    const auto& segment_idx = ctx.collision_pool.particle_indices[rod];

    unsigned sdf_idx = segment_idx[2];
    unsigned rod_idx = ctx.collision_pool.sdf_pool.particles[sdf_idx];
    const auto& rod_params = ctx.collision_pool.sdf_pool.params[sdf_idx];

    // extract current triangle vertex positions
    const Vec3r& v1 = ctx.particles.positions[triangle_idx[0]];
    const Vec3r& v2 = ctx.particles.positions[triangle_idx[1]];
    const Vec3r& v3 = ctx.particles.positions[triangle_idx[2]];

    // extract current rod segment positions
    const Vec3r& s1 = ctx.particles.positions[segment_idx[0]];
    const Vec3r& s2 = ctx.particles.positions[segment_idx[1]];

    Real radius = rod_params.rod.radius;

    Real cp_s;
    Vec3r cp_barys;
    Real dist;
    Math::closestPoint_SegmentTriangle(s1, s2, v1, v2, v3, cp_s, cp_barys, dist);

    if (dist < radius + ctx.params.collision_margin)
    {
        // std::cout << "Rod-triangle collision!" << std::endl;
        Vec3r cp_rod_centerline = s1*(1-cp_s) + s2*cp_s;
        Vec3r cp_tri = v1*cp_barys[0] + v2*cp_barys[1] + v3*cp_barys[2];
        Vec3r diff = cp_tri - cp_rod_centerline;
        Real diff_mag = diff.norm();

        Vec3r normal;
        if (diff_mag > 1e-8)
            normal = diff/diff_mag;
        else
            normal = Vec3r(1,0,0);

        Vec3r cp_rod_global = cp_rod_centerline + normal*radius;
        const Quaternion& q1 = ctx.particles.rotation(segment_idx[0]);
        const Quaternion& q2 = ctx.particles.rotation(segment_idx[1]);
        Quaternion q_mid = Math::Plus_S3(q1, cp_s*Math::Minus_S3(q2, q1));
        Vec3r cp_rod_local = q_mid.inverse() * (radius*normal);

        // std::cout << "Normal: " << normal.transpose() << "  Cp rod local: " << cp_rod_local.transpose() << "  Cp tri: " << cp_tri << std::endl;

        DetectedCollision collision{};
        collision.type = DetectedCollisionType::TriangleRod;
        collision.key = DetectedCollision::generateKey(
            DetectedCollisionType::TriangleRod,
            triangle,
            rod,
            0
        );
        collision.gen1 = ctx.collision_pool.generation[triangle];
        collision.gen2 = ctx.collision_pool.generation[rod];
        collision.normal = normal;
        collision.TriangleRod.tri = Vec3u(triangle_idx[0], triangle_idx[1], triangle_idx[2]);
        collision.TriangleRod.rod = Vec2u(segment_idx[0], segment_idx[1]);
        collision.TriangleRod.barys = cp_barys;
        collision.TriangleRod.s = cp_s;
        collision.TriangleRod.cp_rod_local = cp_rod_local;

        detected_collisions.push_back(std::move(collision));
    }

}

void CollisionDetector::_rodRod(Sim::SimulationContext& ctx, unsigned rod1, unsigned rod2, std::vector<DetectedCollision>& detected_collisions)
{
    /** TODO: (08/24/26) rod-rod collision detection */
}

void CollisionDetector::_sphereSphere(Sim::SimulationContext& ctx, unsigned sphere1, unsigned sphere2, std::vector<DetectedCollision>& detected_collisions)
{
    /** TODO: (07/20/26) sphere-sphere collision detection */
    throw std::runtime_error("Sphere-sphere collision detetction not implemented.");
}

bool CollisionDetector::_triangleSDF_CCD(Sim::SimulationContext& ctx, unsigned triangle, unsigned rb, Real dt, Vec3r& normal, Vec3r& cp_barys, Vec3r& cp_rb_local)
{
    unsigned sdf_idx = ctx.collision_pool.particle_indices[rb][0];      // index in the SDF pool for the rigid body
    unsigned rb_idx = ctx.collision_pool.sdf_pool.particles[sdf_idx];   // particle index of the rigid body COM
    const SDFShapeParams& sdf_params = ctx.collision_pool.sdf_pool.params[sdf_idx];     // SDF parameters

    // extract quantities for rigid body
    const Vec3r& rb_pos = ctx.particles.positions[rb_idx];
    const Quaternion& rb_rot = ctx.particles.rotation(rb_idx);
    const Vec3r& rb_lin_vel = ctx.particles.velocities[rb_idx];
    const Vec3r& rb_ang_vel = ctx.particles.angularVelocity(rb_idx);    // note: body-frame angular velocity

    // evolves rigid body position forward by t
    auto rigid_body_pos = [&](Real t)
    {
        return rb_pos + rb_lin_vel * t;
    };

    // evolves rigid body rotation forward by t
    auto rigid_body_rot = [&](Real t)
    {
        return rb_rot * Math::Exp_s3(rb_ang_vel*t);
    };

    // computes the linear velocity of a point in space if it were attached to the rigid body
    auto lin_vel_on_rigid_body = [&](const Vec3r& xt, Real t)
    {
        const Vec3r rb_ang_vel_global = rigid_body_rot(t) * rb_ang_vel;
        return rb_lin_vel + rb_ang_vel_global.cross(xt - rigid_body_pos(t));
    };

    // computes the relative velocity of a point in space with respect to the rigid body (in the rigid body's local frame)
    auto rel_vel_wrt_rigid_body = [&](const Vec3r& xt, const Vec3r& vt, Real t)
    {
        const Vec3r v_rel_global = vt - lin_vel_on_rigid_body(xt, t);
        return rigid_body_rot(t).inverse() * v_rel_global;
    };

    // transforms a position into the rigid body local frame
    auto local_pos_wrt_rigid_body = [&](const Vec3r& xt, Real t)
    {
        return rigid_body_rot(t).inverse() * (xt - rigid_body_pos(t));
    };

    // extract quantities for triangle
    const auto& tri_indices = ctx.collision_pool.particle_indices[triangle];
    const Vec3r& tri1 = ctx.particles.positions[tri_indices[0]];
    const Vec3r& tri2 = ctx.particles.positions[tri_indices[1]];
    const Vec3r& tri3 = ctx.particles.positions[tri_indices[2]];
    const Vec3r& tri_vel1 = ctx.particles.velocities[tri_indices[0]];
    const Vec3r& tri_vel2 = ctx.particles.velocities[tri_indices[1]];
    const Vec3r& tri_vel3 = ctx.particles.velocities[tri_indices[2]];

    /** Find initial iterate */
    // helper for evaluating which vertex of the triangle to start at
    // use the vertex that minimizes vi^T * gradSDF(si) where si are the triangle vertices and vi the corresponding linear velocities
    auto evaluate_initial_iterate = [&](const Vec3r& si, const Vec3r& vi)
    {
        const Vec3r rel_vel = rel_vel_wrt_rigid_body(si, vi, 0);
        const Vec3r rel_pos = local_pos_wrt_rigid_body(si, 0);

        // vi^T * gradSDF(si) (everything expressed in SDF local frame)
        return rel_vel.dot(SDF::gradient(sdf_params, rel_pos));
    };

    // evaluate each vertex
    Real eval1 = evaluate_initial_iterate(tri1, tri_vel1);
    Real eval2 = evaluate_initial_iterate(tri2, tri_vel2);
    Real eval3 = evaluate_initial_iterate(tri3, tri_vel3);

    // optimization variables
    Real t_start = 0, t_end = dt;  // time interval
    Real t = 0; // time
    Real u,v,w; // barycentric coordinates on triangle

    // vertex 1 is the initial iterate
    if (eval1 < eval2 && eval1 < eval3)         { u = 1; v = 0; w = 0; }
    // vertex 2 is the initial iterate
    else if (eval2 < eval1 && eval2 < eval3)    { u = 0; v = 1; w = 0; }
    // vertex 3 is the initial iterate
    else                                        { u = 0; v = 0; w = 1; }


    /** Spatio-temporal optimization */
    // helper that interpolates the barycentric position in time
    // note: this is still the global position - will need to convert to local SDF frame
    auto barycentric_interpolate = [&](Real u, Real v, Real w, Real t)
    {
        return u*(tri1 + tri_vel1*t) + v*(tri2 + tri_vel2*t) + w*(tri3 + tri_vel3*t);
    };
    // helper for velocity interpolation of point on triangle (global frame)
    auto barycentric_interpolate_velocity = [&](Real u, Real v, Real w)
    {
        return u*tri_vel1 + v*tri_vel2 + w*tri_vel3;
    };
    // evaluate the SDF at a given time
    auto signed_distance_at_time = [&](Real t)
    {
        Vec3r x_global = barycentric_interpolate(u, v, w, t);
        Vec3r x = local_pos_wrt_rigid_body(x_global, t);
        return SDF::evaluate(sdf_params, x);
    };
    // evaluate the abs of the SDF at a given time
    auto unsigned_distance_at_time = [&](Real t)
    {
        return std::abs(signed_distance_at_time(t));
    };

    // evaluate SDF at a given point - assumes positions are in rigid body local frame
    auto signed_distance_at_point = [&](const Vec3r& x_local)
    {
        return SDF::evaluate(sdf_params, x_local);
    };

    constexpr unsigned MAX_ITER = 16;
    constexpr Real TOL = 1e-6;
    Real sdf_xti;   // latest SDF evaluation
    Vec3r grad_xti; // latest SDF gradient
    Vec3r x_ti;     // latest contact point on triangle
    for (unsigned iter = 0; iter < MAX_ITER; iter++)
    {
        /** Solve temporal subproblem */
        Real t_old = t;
        Vec3r x_ti_global = barycentric_interpolate(u, v, w, t);
        Vec3r v_ti_global = barycentric_interpolate_velocity(u, v, w);
        x_ti = local_pos_wrt_rigid_body(x_ti_global, t);
        Vec3r v_ti = rel_vel_wrt_rigid_body(x_ti_global, v_ti_global, t);

        sdf_xti = SDF::evaluate(sdf_params, x_ti);
        if (sdf_xti <= 0)
        {
            // update the interval
            t_end = std::min(t, t_end);
            t = Algorithm::goldenSectionSearch(t_start, t, unsigned_distance_at_time);
        }
        else
        {
            Real v_dot_grad = v_ti.dot(SDF::gradient(sdf_params, x_ti));
            Real d = v_dot_grad > 0 ? -1 : 1;
            if (d < 0)
                t = Algorithm::goldenSectionSearch(t_start, t, signed_distance_at_time);
            else
                t = Algorithm::goldenSectionSearch(t, t_end, signed_distance_at_time);
        }

        /** Solve spatial subproblem */
        x_ti_global = barycentric_interpolate(u, v, w, t);
        v_ti_global = barycentric_interpolate_velocity(u, v, w);
        x_ti = local_pos_wrt_rigid_body(x_ti_global, t);
        v_ti = rel_vel_wrt_rigid_body(x_ti_global, v_ti_global, t);

        sdf_xti = SDF::evaluate(sdf_params, x_ti);
        grad_xti = SDF::gradient(sdf_params, x_ti);
        if (sdf_xti <= 0)
            t_end = std::min(t, t_end);
        
        // find support vertex
        Vec3r tri1_localt = local_pos_wrt_rigid_body(tri1 + tri_vel1*t, t);
        Vec3r tri2_localt = local_pos_wrt_rigid_body(tri2 + tri_vel2*t, t);
        Vec3r tri3_localt = local_pos_wrt_rigid_body(tri3 + tri_vel3*t, t);
        Real s1 = tri1_localt.dot(grad_xti);
        Real s2 = tri2_localt.dot(grad_xti);
        Real s3 = tri3_localt.dot(grad_xti);
        Vec3r s;
        if (s1 < s2 && s1 < s3)         s = tri1_localt;
        else if (s2 < s1 && s2 < s3)    s = tri2_localt;
        else                            s = tri3_localt;

        Vec3r x_new = Algorithm::goldenSectionSearch(x_ti, s, signed_distance_at_point);
        Vec3r barys = Math::barycentricCoordinates(x_new, tri1_localt, tri2_localt, tri3_localt);
        u = barys[0]; v = barys[1]; w = barys[2];

        // see if we can exit early
        if (std::abs(t_old - t) < TOL && (x_new - x_ti).norm() < TOL)
            break;
    }

    cp_barys = Vec3r(u, v, w);
    cp_rb_local = x_ti - sdf_xti*grad_xti;
    normal = rigid_body_rot(t) * grad_xti;

    // std::cout << "Best t: " << t << std::endl;
    // std::cout << "SDF @ t: " << signed_distance_at_time(t) << std::endl;

    return sdf_xti < ctx.params.collision_margin;
    

}

} // namespace Collision