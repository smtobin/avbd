#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "collision/AABB.hpp"

namespace Collision
{

enum class PrimitiveType : uint8_t {
    RodSegment,
    Triangle,
    RigidSphere,
    RigidBox,
    RigidSDF
};

/** Storage of different collision primitives. Used in collision detection and generating the LBVH. */
/** TODO: Think about integration of oriented particles (07/16/26) */
struct CollisionPrimitivePool : TombstonePool
{
    /** Collision primitive storage */
    std::vector<PrimitiveType> type;
    std::vector<std::array<unsigned, 3>> particle_indices;
    std::vector<uint8_t> num_particles;
    std::vector<unsigned> object_id;

    /** Collision geometry, recomputed every frame */
    std::vector<AABB> aabb;
    std::vector<Vec3r> centroid;
    std::vector<uint64_t> morton_code;

    /** BVH bookkeeping */
    std::vector<unsigned> sorted_order;

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit CollisionPrimitivePool(unsigned capacity)
        : TombstonePool(capacity)
        , type(capacity)
        , particle_indices(capacity)
        , num_particles(capacity)
        , object_id(capacity)
        , aabb(capacity)
        , centroid(capacity)
        , morton_code(capacity)
        , sorted_order(capacity)
    {

    }

    /** Add objects
     * 
     * TODO: (07/17/26) Make this store object id and a pointer to update with changing topology
     * Not really necessary right now
     */
    void addObject(const SimObject::TetMeshObject& mesh_obj);
    void addObject(const SimObject::RigidSphere& sphere);
};

} // namespace Collision