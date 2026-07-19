#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "collision/AABB.hpp"
#include "collision/SDFPrimitivePool.hpp"

namespace Collision
{

enum class PrimitiveType : uint8_t {
    RodSegment,
    Triangle,
    RigidSDF
};


/** Storage of different collision primitives. Used in collision detection and generating the LBVH. */
/** TODO: Think about integration of oriented particles (07/16/26) */
struct CollisionPrimitivePool : TombstonePool
{
    /** Collision primitive storage */
    // type of collision primitive
    std::vector<PrimitiveType> type;
    // storage of up to 3 particle indices
    // if the PrimitveType == RigidSDF, then particle_indices[0] is the index in the SDFPrimitivePool
    std::vector<std::array<unsigned, 3>> particle_indices;
    std::vector<uint8_t> num_particles;     // number of particles for each primitive (not sure if this is really necessary)
    std::vector<unsigned> object_id;        /** TODO: (07/19/26) do something with object ID */

    SDFPrimitivePool sdf_pool;  // storage of extra data for SDFs

    /** Collision geometry, recomputed every frame */
    std::vector<AABB> aabb;
    std::vector<Vec3r> centroid;
    std::vector<uint64_t> morton_code;

    /** BVH bookkeeping */
    std::vector<unsigned> sorted_order;

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool for total primitives
     * @param sdf_capacity : the capacity of the memory pool for objects described by an SDF
     */
    explicit CollisionPrimitivePool(unsigned capacity, unsigned sdf_capacity)
        : TombstonePool(capacity)
        , type(capacity)
        , particle_indices(capacity)
        , num_particles(capacity)
        , object_id(capacity)
        , sdf_pool(sdf_capacity)
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