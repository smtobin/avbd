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
    // build BVH
    LBVHBuilder::buildBVH(ctx.particles, ctx.collision_pool, ctx.lbvh);

    // traverse BVH for potential collisions
    LBVHTraversal::traverseSelfIterative(ctx.lbvh, ctx.lbvh.root, _potential_collisions);

    // narrow-phase collision detection
    _narrowPhaseCollisionDetection(ctx);

    // handle detected collisions
    // sort the current detected collisions by key
    // merge lists, and create/destroy collision constraints accordingly
    _handleDetectedCollisions(ctx);
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

    // check distance between closest point and sphere center
    Vec3r diff = (tri_cp - p);
    Real dist = diff.norm();
    if (dist < sphere_params.sphere.radius)
    {
        // compute collision normal
        Vec3r normal;
        if (dist > Real(1e-6))
            normal = diff/dist;
        else
            normal = Vec3r(1,0,0);  // arbitrary direction

        // barycentric coordinates of triangle collision point
        Vec3r barys = Math::barycentricCoordinates(tri_cp, v1, v2, v3);

        // sphere conact point in local frame
        Vec3r sphere_cp_local = ctx.particles.rotation(sphere_idx).inverse() * normal*sphere_params.sphere.radius;

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
        collision.TriangleRigid.barys = barys;
        collision.TriangleRigid.rb = sphere_idx;
        collision.TriangleRigid.cp_rb_local = sphere_cp_local;

        _cur_detected_collisions.push_back(std::move(collision));
    }

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

} // namespace Collision