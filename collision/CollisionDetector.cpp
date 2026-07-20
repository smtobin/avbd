#include "collision/CollisionDetector.hpp"

#include "simulation/SimulationContext.hpp"

#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"

#include "common/Math.hpp"

namespace Collision
{

void CollisionDetector::detectCollisionsAndRecolor(Sim::SimulationContext& ctx)
{
    if (!_collision_table_initialized)
        _initCollisionTable();

    // build BVH
    LBVHBuilder::buildBVH(ctx.particles, ctx.oriented_particles, ctx.collision_pool, ctx.lbvh);

    // traverse BVH for potential collisions
    std::vector<std::pair<unsigned, unsigned>> potential_collisions;
    LBVHTraversal::traverseSelfIterative(ctx.lbvh, ctx.lbvh.root, potential_collisions);

    // narrow-phase collision detection
    _narrowPhaseCollisionDetection(potential_collisions);

}

void CollisionDetector::_narrowPhaseCollisionDetection(Sim::SimulationContext& ctx, const std::vector<std::pair<unsigned, unsigned>>& potential_collisions)
{
    /** TODO: (07/19/26) Parallelize this? */
    for (const auto& potential_collision : potential_collisions)
    {
        /** TODO: (07/19/26) Need to figure this out and get it straight and document it:
         * Potential collisions = pairs of (leaf) nodes in the BVH. Leaf node index <=> primitive index in collision pool?
         * Or is it leaf node index <=> sorted order index ==> primitive index?
         */
        unsigned a = potential_collision.first;
        unsigned b = potential_collision.second;

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

            case default:
            {
                throw std::runtime_error("Unsupported collision!");
            }
            
        }
    }
}

void CollisionDetector::_triangleSphere(Sim::SimulationContext& ctx, unsigned triangle, unsigned sphere)
{
    const auto& triangle_idx = ctx.collision_pool.particle_indices[triangle];

    unsigned sdf_idx = ctx.collision_pool.particle_indices[sphere][0];
    unsigned sphere_idx = ctx.collision_pool.sdf_pool.particles[sdf_idx];
    const auto& sphere_params = ctx.collision_pool.sdf_pool.params[sdf_idx];

    // extract current triangle vertex positions
    const Vec3r& v1 = ctx.particles.positions[triangle_idx[0]];
    const Vec3r& v2 = ctx.particles.positions[triangle_idx[1]];
    const Vec3r& v3 = ctx.particles.positions[triangle_idx[2]];

    // extract current sphere center
    const Vec3r& p = ctx.oriented_particles.positions[sphere_idx];

    // closest point on triangle to sphere center
    Vec3r tri_cp = Math::closestPointOnTriangle(p, v1, v2, v3);

    // check distance between closest point and sphere center
    Real dist = (tri_cp - p).norm();
    if (dist < sphere_params.sphere.radius)
    {
        std::cout << "Sphere-triangle collision!" << std::endl;
    }

}

void CollisionDetector::_triangleTriangle(Sim::SimulationContext& ctx, unsigned triangle1, unsigned triangle2)
{
    
}

void CollisionDetector::_sphereSphere(Sim::SimulationContext& ctx, unsigned sphere1, unsigned sphere2)
{

}

} // namespace Collision