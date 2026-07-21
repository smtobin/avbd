#include "collision/CollisionDetector.hpp"

#include "simulation/SimulationContext.hpp"

#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"

#include "common/Math.hpp"

namespace Collision
{

void CollisionDetector::detectCollisionsAndRecolor(Sim::SimulationContext& ctx)
{
    // build BVH
    LBVHBuilder::buildBVH(ctx.particles, ctx.oriented_particles, ctx.collision_pool, ctx.lbvh);

    // traverse BVH for potential collisions
    std::vector<std::pair<unsigned, unsigned>> potential_collisions;
    LBVHTraversal::traverseSelfIterative(ctx.lbvh, ctx.lbvh.root, potential_collisions);

    // narrow-phase collision detection
    _narrowPhaseCollisionDetection(ctx, potential_collisions);

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

void CollisionDetector::_narrowPhaseCollisionDetection(Sim::SimulationContext& ctx, const std::vector<std::pair<unsigned, unsigned>>& potential_collisions)
{
    /** TODO: (07/19/26) Parallelize this? */
    for (const auto& potential_collision : potential_collisions)
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
    const Vec3r& p = ctx.oriented_particles.positions[sphere_idx];

    // closest point on triangle to sphere center
    Vec3r tri_cp = Math::closestPoint_PointTriangle(p, v1, v2, v3);

    // check distance between closest point and sphere center
    Real dist = (tri_cp - p).norm();
    if (dist < sphere_params.sphere.radius)
    {
        std::cout << "Sphere-triangle collision!" << std::endl;
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