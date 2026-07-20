#include "collision/CollisionDetector.hpp"

#include "simulation/SimulationContext.hpp"

#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"

namespace Collision
{

// initialize static values
bool CollisionDetector::_collision_table_initialized = false;
CollisionDetector::CollisionFunc CollisionDetector::_collision_table[static_cast<int>(ColliderType::COUNT)][static_cast<int>(ColliderType::COUNT)] = {};

// [](void* a, void* b) {
    // CollisionScene::_checkCollision(static_cast<SimObject::XPBDRigidSphere*>(a), static_cast<SimObject::XPBDRigidSphere*>(b));
// };
void CollisionDetector::_initCollisionTable()
{
    if (_collision_table_initialized)
        return;
    
    // first type is a triangle in a mesh
    _collision_table[static_cast<unsigned>(PrimitiveType::Triangle)][static_cast<unsigned>(PrimitiveType::Triangle)] = CollisionDetector::_triangleTriangleCollision;

    

    _collision_table_initialized = true;
}

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
    /** TODO: (07/19/26) Parallelize this? */
    for (const auto& potential_collision : potential_collisions)
    {
        /** TODO: (07/19/26) Need to figure this out and get it straight and document it:
         * Potential collisions = pairs of (leaf) nodes in the BVH. Leaf node index <=> primitive index in collision pool?
         * Or is it leaf node index <=> sorted order index ==> primitive index?
         */
        PrimitiveType type1 = ctx.collision_pool.type[potential_collision.first];
        PrimitiveType type2 = ctx.collision_pool.type[potential_collision.second];


    }

}

} // namespace Collision