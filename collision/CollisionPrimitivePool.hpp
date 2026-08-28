#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "common/ParticlePool.hpp"
#include "collision/AABB.hpp"
#include "collision/SDFPrimitivePool.hpp"
#include "simulation/SimulationParams.hpp"

namespace Collision
{

/** The type of primitive
 * 
 * This is used for storage in the CollisionPrimitivePool.
 * 
 * RodSegment: a segment of a rod, generally one element (not implemented)
 * Triangle: a triangular face in a mesh - defined by the 3 vertices
 * RigidSDF: a rigid object represented by its SDF. See SDFType enum for supported SDF types.
 */
enum class PrimitiveType : uint8_t {
    RodSegment,
    Triangle,
    RigidSDF
};

/** The type of collision geometry
 * 
 * This is used for dispatching the proper narrow-phase collision subroutine.
 * Basically just makes a flat list of supported geometry in the sim, including the types of rigid SDFs supported.
 */
enum class CollisionGeometryType : uint8_t {
    Triangle=0,
    RodSegment,

    // rigid SDF types
    Sphere,
    Box,
    Capsule,

    // total number of collision geometries
    Count
};


/** Storage of different collision primitives. Used in collision detection and generating the LBVH. */
struct CollisionPrimitivePool : TombstonePool
{
    /** Collision primitive storage */
    // type of collision primitive
    std::vector<PrimitiveType> type;
    // storage of up to 3 particle indices
    // if the PrimitveType == RigidSDF, then particle_indices[0] is the index in the SDFPrimitivePool
    // if the PrimitiveType == Rod, then particle_indices[2] is the the index in the SDFPrimitivePool
    std::vector<std::array<unsigned, 3>> particle_indices;
    std::vector<uint8_t> num_particles;     // number of particles for each primitive (not sure if this is really necessary)
    std::vector<unsigned> object_id;

    /** TODO: (08/24/26)
     * Rename SDFPrimitivePool to reflect general storage of auxiliary data for collision primitives
     * e.g. "PrimitiveAuxData"
     */
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

    /** Get the CollisionGeometryType for a primitive in the pool */
    inline CollisionGeometryType getCollisionGeometryType(unsigned p_idx)
    {
        switch(type[p_idx])
        {
            case PrimitiveType::RodSegment:
                return CollisionGeometryType::RodSegment;
            case PrimitiveType::Triangle:
                return CollisionGeometryType::Triangle;
            case PrimitiveType::RigidSDF:
                switch(sdf_pool.params[particle_indices[p_idx][0]].type)
                {
                    case SDFType::Sphere:
                        return CollisionGeometryType::Sphere;
                    case SDFType::Box:
                        return CollisionGeometryType::Box;
                    case SDFType::Capsule:
                        return CollisionGeometryType::Capsule;
                    default:
                        throw std::runtime_error("Unrecognized SDF type");
                }
            default:
                throw std::runtime_error("Unrecognized collision geometry type");
        }
    }

    /** AABB for an object */
    inline AABB globalBounds(unsigned p_idx, const ParticlePool& particle_pool) const
    {
        switch(type[p_idx])
        {
            case PrimitiveType::Triangle:
            {
                AABB box = AABB::empty();
                for (unsigned k = 0; k < 3; k++)
                {
                    // expand box based on each triangle vertex
                    box.expand(particle_pool.positions[particle_indices[p_idx][k]]);
                }
                return box;
            }
            case PrimitiveType::RigidSDF:
            {
                // index of the SDF params in the SDF pool
                unsigned sdf_idx = particle_indices[p_idx][0];
                AABB local = SDF::localBounds(sdf_pool.params[sdf_idx]);
                
                // transform local AABB to global
                Vec3r center = 0.5 * (local.min + local.max);
                Vec3r extent = 0.5 * (local.max - local.min);

                // index of the oriented particle in the oriented particle pool
                unsigned op_idx = sdf_pool.particles[sdf_idx];
                const Quaternion& rotation = particle_pool.rotation(op_idx);
                Vec3r world_center = rotation * center + particle_pool.positions[op_idx];
                
                Mat3r abs_R = rotation.toRotationMatrix().cwiseAbs();
                Vec3r world_extent = abs_R * extent;

                return { world_center - world_extent,  world_center + world_extent };
            }
            case PrimitiveType::RodSegment:
            {
                const Vec3r& p1 = particle_pool.positions[particle_indices[p_idx][0]];
                const Vec3r& p2 = particle_pool.positions[particle_indices[p_idx][1]];
                Real radius = sdf_pool.params[particle_indices[p_idx][2]].rod.radius;
                Vec3r bbox_min = p1.cwiseMin(p2) - Vec3r::Constant(radius);
                Vec3r bbox_max = p1.cwiseMax(p2) + Vec3r::Constant(radius);
                
                return { bbox_min, bbox_max };
            }
            default:
                throw std::runtime_error("Unrecognized primitive type");
        }
    }

    /** AABB for an object, expanded by its current velocity. Useful for predictive collision detection */
    inline AABB speculativeGlobalBounds(unsigned p_idx, const ParticlePool& particle_pool, Real dt)
    {
        switch(type[p_idx])
        {
            case PrimitiveType::Triangle:
            {
                AABB box = AABB::empty();
                for (unsigned k = 0; k < 3; k++)
                {
                    // expand box based on each triangle vertex
                    box.expand(particle_pool.positions[particle_indices[p_idx][k]]);
                    box.expand(particle_pool.positions[particle_indices[p_idx][k]] + dt * particle_pool.velocities[particle_indices[p_idx][k]]);
                }
                return box;
            }
            case PrimitiveType::RigidSDF:
            {
                // index of the SDF params in the SDF pool
                unsigned sdf_idx = particle_indices[p_idx][0];
                AABB local = SDF::localBounds(sdf_pool.params[sdf_idx]);
                
                // transform local AABB to global
                Vec3r center = 0.5 * (local.min + local.max);
                Vec3r halfextent = 0.5 * (local.max - local.min);

                // index of the oriented particle in the oriented particle pool
                unsigned op_idx = sdf_pool.particles[sdf_idx];

                // compute current AABB
                const Quaternion& rotation = particle_pool.rotation(op_idx);
                Vec3r world_center = rotation * center + particle_pool.positions[op_idx];
                Mat3r abs_R = rotation.toRotationMatrix().cwiseAbs();
                Vec3r world_extent = abs_R * halfextent;
                AABB box = { world_center - world_extent, world_center + world_extent };

                // compute AABB at predicted position and rotation given the current linear and angular velocities
                const Quaternion new_rotation = rotation; /** TODO: (07/22/26) Predict new orientation using exp map or linearized update */
                Vec3r new_world_center = new_rotation * center + particle_pool.positions[op_idx] + dt * particle_pool.velocities[op_idx];
                Mat3r new_abs_R = new_rotation.toRotationMatrix().cwiseAbs();
                Vec3r new_world_extent = new_abs_R * halfextent;
                AABB pred_box = { new_world_center - new_world_extent, new_world_center + new_world_extent };

                // merge the boxes
                box.expand(pred_box);

                return box;
            }
            case PrimitiveType::RodSegment:
            {
                const Vec3r& p1 = particle_pool.positions[particle_indices[p_idx][0]];
                const Vec3r& p2 = particle_pool.positions[particle_indices[p_idx][1]];
                Real radius = sdf_pool.params[particle_indices[p_idx][2]].rod.radius;
                Vec3r bbox_min = p1.cwiseMin(p2) - Vec3r::Constant(radius);
                Vec3r bbox_max = p1.cwiseMax(p2) + Vec3r::Constant(radius);
                AABB bbox = { bbox_min, bbox_max };

                Vec3r next_p1 = p1 + dt * particle_pool.velocities[particle_indices[p_idx][0]];
                Vec3r next_p2 = p2 + dt * particle_pool.velocities[particle_indices[p_idx][1]];
                Vec3r next_bbox_min = next_p1.cwiseMin(next_p2) - Vec3r::Constant(radius);
                Vec3r next_bbox_max = next_p1.cwiseMax(next_p2) + Vec3r::Constant(radius);
                AABB next_bbox = { next_bbox_min, next_bbox_max };

                // merge the boxes
                bbox.expand(next_bbox);
                return bbox;
            }
            default:
                throw std::runtime_error("Unrecognized primitive type");
        }
    }

    /** Add objects
     * 
     * TODO: (07/17/26) Make this store object id and a pointer to update with changing topology
     * Not really necessary right now
     */
    void addObject(const SimObject::TetMeshObject& mesh_obj);
    unsigned addObject(const SimObject::RigidSphere& sphere);

    /** @returns the index in the SDF pool */
    unsigned addObject(const SimObject::Rod& rod);
};

} // namespace Collision