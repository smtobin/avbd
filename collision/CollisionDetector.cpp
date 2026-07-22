#include "collision/CollisionDetector.hpp"

#include "simulation/SimulationContext.hpp"

#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"

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
    LBVHBuilder::buildBVH(ctx.particles, ctx.collision_pool, ctx.lbvh, ctx.params.dt);

    // traverse BVH for potential collisions
    LBVHTraversal::traverseSelfIterative(ctx.lbvh, ctx.lbvh.root, _potential_collisions);

    // narrow-phase collision detection
    _narrowPhaseCollisionDetection(ctx);

    // handle detected collisions
    // sort the current detected collisions by key
    // merge lists, and create/destroy collision constraints accordingly
    _handleDetectedCollisions(ctx);

    // rebuild adjacency
    /** TODO: (07/22/26) Do this incrementally? Handle collision constraints separately? Anything to not recompute adjacency from scratch each time. */
    ctx.adjacency.buildAdjacency(ctx.particles, ctx.energies);
    ctx.coloring.buildColorList(ctx.adjacency, ctx.particles.totalSize());
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

void CollisionDetector::_narrowPhaseCollisionDetection(Sim::SimulationContext& ctx)
{
    /** TODO: (07/19/26) Parallelize this? */
    for (const auto& potential_collision : _potential_collisions)
    {
        // potential collisions are pairs of LBVH leaf node indices, which correspond to the collision pool's sorted order
        // extract the primitive indices by indexing in the sorted order array
        unsigned a = ctx.collision_pool.sorted_order[potential_collision.first];
        unsigned b = ctx.collision_pool.sorted_order[potential_collision.second];

        if (_shouldSkip(ctx.collision_pool, a, b))
            continue;

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
                _triangleTriangle(ctx, a, b);
                break;
            }
            case _makeCollisionKey(CollisionGeometryType::Triangle, CollisionGeometryType::Sphere):
            {
                _triangleSphere(ctx, a, b);
                break;
            }
            case _makeCollisionKey(CollisionGeometryType::Sphere, CollisionGeometryType::Sphere):
            {
                _sphereSphere(ctx, a, b);
                break;
            }

            default:
            {
                throw std::runtime_error("Unsupported collision!");
            }
            
        }
    }
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
            unsigned slot = ctx.energies.triangle_rigid_collision.addEnergy(
                collision.TriangleRigid.tri[0],
                collision.TriangleRigid.tri[1],
                collision.TriangleRigid.tri[2],
                collision.TriangleRigid.rb,
                nullptr,
                collision.normal,
                collision.TriangleRigid.barys,
                collision.TriangleRigid.cp_rb_local
            );
            collision.e_idx = slot;
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
            info.normal = collision.normal;
            info.cp_rb_local = collision.TriangleRigid.cp_rb_local;
            info.barys = collision.TriangleRigid.barys;
            break;
        }
        default:
        {
            throw std::runtime_error("CollisionDetector::_updateCollision - unsupported DetectedCollisionType!");
            break;
        }
    }
}


void CollisionDetector::_triangleSphere(Sim::SimulationContext& ctx, unsigned triangle, unsigned sphere)
{
    std::cout << "Testing sphere-triangle collision..." << std::endl;
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

        _cur_detected_collisions.push_back(std::move(collision));
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

void CollisionDetector::_triangleTriangle(Sim::SimulationContext& ctx, unsigned triangle1, unsigned triangle2)
{
    /** TODO: (07/20/26) triangle-triangle collision detection */
    // std::cout << "Testing triangle-triangle collision..." << std::endl;
    // throw std::runtime_error("Triangle-triangle collision detection not implemented.");
}

void CollisionDetector::_sphereSphere(Sim::SimulationContext& ctx, unsigned sphere1, unsigned sphere2)
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
        return rb_rot * Math::QuaternionExp_so3(rb_ang_vel*t);
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

    std::cout << "Best t: " << t << std::endl;
    std::cout << "SDF @ t: " << signed_distance_at_time(t) << std::endl;

    return sdf_xti < ctx.params.collision_margin;
    

}

} // namespace Collision