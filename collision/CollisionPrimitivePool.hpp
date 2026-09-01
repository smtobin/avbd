#pragma once

#include "common/common.hpp"
#include "common/Math.hpp"
#include "common/TombstonePool.hpp"
#include "common/ParticlePool.hpp"
#include "collision/CollisionGeometryType.hpp"
#include "collision/CollisionShapeParams.hpp"
#include "collision/SDF.hpp"
#include "collision/AABB.hpp"
#include "simulation/SimulationParams.hpp"

namespace Collision
{

/** Storage of extra data for analytic collision shapes. */
struct CollisionShapeParamsPool : TombstonePool
{
    std::vector<CollisionShapeParams> shape_params;

    explicit CollisionShapeParamsPool(unsigned capacity)
        : TombstonePool(capacity)
        , shape_params(capacity)
    {}

    /** Create shape params for an object.
     * @returns the index
     */
    unsigned addObject(const SimObject::RigidSphere& sphere);
    unsigned addObject(const SimObject::Rod& rod);

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Storage of different collision primitives. Used in collision detection and generating the LBVH. */
struct CollisionPrimitivePool : TombstonePool
{
    /** Collision primitive storage */
    // type of collision primitive
    std::vector<CollisionGeometryType> type;
    // storage of up to 3 particle indices
    // if the type is a rigid body (sphere, box, etc.), then particle_indices[0] is the COM index and particle_indices[1] is index into the shape params vector
    // if the type == RodSegment, then particle_indices[2] is the the index in the shape params vector
    // this indexing is assisted by the getShapeParamIndex() helper function
    std::vector<std::array<unsigned, 3>> particle_indices;
    std::vector<uint8_t> num_particles;     // number of particles for each primitive (not sure if this is really necessary)
    std::vector<unsigned> object_id;

    // required data for resolving collisions with analytic shapes (e.g. spheres, boxes)
    // one entry per analytic object (not per-primitive)
    CollisionShapeParamsPool shape_params_pool;

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
        , shape_params_pool(sdf_capacity)
        , aabb(capacity)
        , centroid(capacity)
        , morton_code(capacity)
        , sorted_order(capacity)
    {

    }

    /** Helper for returning the index in the shape params vector for a given primitive.
     * Assumes that the given primitive actually has extra data associated with it.
     * If it does not (e.g., it is a triangle primitive), this function throws an exception.
     */
    inline unsigned shapeParamsIndex(unsigned p_idx) const
    {
        if (isAnalyticShape(type[p_idx]))
        {
            // for analytic shapes, [0] is the COM particle and [1] is the index in the shape params vector
            return particle_indices[p_idx][1];
        }
        else if (type[p_idx] == CollisionGeometryType::RodSegment)
        {
            // for rod segments, [0] and [1] are the rod nodes, and [2] is the index in the shape params vector
            return particle_indices[p_idx][2];
        }
        else
        {
            // unsupported type that does not actually have extra data
            throw std::runtime_error("CollisionGeometryType does not have extra data associated with it.");
        }
    }

    /** AABB for an object */
    inline AABB globalBounds(unsigned p_idx, const ParticlePool& particle_pool) const
    {
        if (type[p_idx] == CollisionGeometryType::Triangle)
        {
            AABB box = AABB::empty();
            for (unsigned k = 0; k < 3; k++)
            {
                // expand box based on each triangle vertex
                box.expand(particle_pool.positions[particle_indices[p_idx][k]]);
            }
            return box;
        }
        else if (type[p_idx] == CollisionGeometryType::RodSegment)
        {
            // get rod node positions
            const Vec3r& p1 = particle_pool.positions[particle_indices[p_idx][0]];
            const Vec3r& p2 = particle_pool.positions[particle_indices[p_idx][1]];

            // get rod radius
            unsigned params_idx = shapeParamsIndex(p_idx);
            Real radius = shape_params_pool.shape_params[params_idx].rod.radius;

            // construct bounding box
            Vec3r bbox_min = p1.cwiseMin(p2) - Vec3r::Constant(radius);
            Vec3r bbox_max = p1.cwiseMax(p2) + Vec3r::Constant(radius);
            
            return { bbox_min, bbox_max };
        }
        else if (isAnalyticShape(type[p_idx]))
        {
            // index of the SDF params in the SDF pool
            unsigned params_idx = shapeParamsIndex(p_idx);
            AABB local = SDF::localBounds(shape_params_pool.shape_params[params_idx]);
            
            // transform local AABB to global
            Vec3r center = 0.5 * (local.min + local.max);
            Vec3r extent = 0.5 * (local.max - local.min);

            // index of the oriented particle in the oriented particle pool
            unsigned op_idx = particle_indices[p_idx][0];
            const Quaternion& rotation = particle_pool.rotation(op_idx);
            Vec3r world_center = rotation * center + particle_pool.positions[op_idx];
            
            Mat3r abs_R = rotation.toRotationMatrix().cwiseAbs();
            Vec3r world_extent = abs_R * extent;

            return { world_center - world_extent,  world_center + world_extent };
        }
        else
        {
            throw std::runtime_error("globalBounds not implemented for this CollisionGeometryType.");
        }
    }

    /** AABB for an object, expanded by its current velocity. Useful for predictive collision detection */
    inline AABB speculativeGlobalBounds(unsigned p_idx, const ParticlePool& particle_pool, Real dt)
    {
        if (type[p_idx] == CollisionGeometryType::Triangle)
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
        else if (type[p_idx] == CollisionGeometryType::RodSegment)
        {
            // get rod node positions
            const Vec3r& p1 = particle_pool.positions[particle_indices[p_idx][0]];
            const Vec3r& p2 = particle_pool.positions[particle_indices[p_idx][1]];

            // get radius from shape params
            unsigned params_idx = shapeParamsIndex(p_idx);
            Real radius = shape_params_pool.shape_params[params_idx].rod.radius;

            // construct current bounding box
            Vec3r bbox_min = p1.cwiseMin(p2) - Vec3r::Constant(radius);
            Vec3r bbox_max = p1.cwiseMax(p2) + Vec3r::Constant(radius);
            AABB bbox = { bbox_min, bbox_max };

            // construct extrapolated bounding box at next frame based on current velocities
            Vec3r next_p1 = p1 + dt * particle_pool.velocities[particle_indices[p_idx][0]];
            Vec3r next_p2 = p2 + dt * particle_pool.velocities[particle_indices[p_idx][1]];
            Vec3r next_bbox_min = next_p1.cwiseMin(next_p2) - Vec3r::Constant(radius);
            Vec3r next_bbox_max = next_p1.cwiseMax(next_p2) + Vec3r::Constant(radius);
            AABB next_bbox = { next_bbox_min, next_bbox_max };

            // merge the boxes
            bbox.expand(next_bbox);
            return bbox;
        }
        else if (isAnalyticShape(type[p_idx]))
        {
            // index of the SDF params in the SDF pool
            unsigned params_idx = shapeParamsIndex(p_idx);
            AABB local = SDF::localBounds(shape_params_pool.shape_params[params_idx]);
            
            // transform local AABB to global
            Vec3r center = 0.5 * (local.min + local.max);
            Vec3r halfextent = 0.5 * (local.max - local.min);

            // index of the oriented particle in the oriented particle pool
            unsigned op_idx = particle_indices[p_idx][0];

            // compute current AABB
            const Quaternion& rotation = particle_pool.rotation(op_idx);
            Vec3r world_center = rotation * center + particle_pool.positions[op_idx];
            Mat3r abs_R = rotation.toRotationMatrix().cwiseAbs();
            Vec3r world_extent = abs_R * halfextent;
            AABB box = { world_center - world_extent, world_center + world_extent };

            // compute AABB at predicted position and rotation given the current linear and angular velocities
            const Vec3r& ang_vel = particle_pool.angularVelocity(op_idx);
            const Quaternion new_rotation = Math::Plus_S3(rotation, dt*ang_vel);
            Vec3r new_world_center = new_rotation * center + particle_pool.positions[op_idx] + dt * particle_pool.velocities[op_idx];
            Mat3r new_abs_R = new_rotation.toRotationMatrix().cwiseAbs();
            Vec3r new_world_extent = new_abs_R * halfextent;
            AABB pred_box = { new_world_center - new_world_extent, new_world_center + new_world_extent };

            // merge the boxes
            box.expand(pred_box);

            return box;
        }
        else
        {
            throw std::runtime_error("speculativeGlobalBounds not implemented for this CollisionGeometryType.");
        }
    }

    /** Add objects */
    void addObject(const SimObject::TetMeshObject& mesh_obj);
    unsigned addObject(const SimObject::RigidSphere& sphere);
    unsigned addObject(const SimObject::Rod& rod);
};

} // namespace Collision