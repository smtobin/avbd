#pragma once

#include "common/common.hpp"
#include "collision/CollisionPrimitivePool.hpp"
#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"
#include "collision/DetectedCollision.hpp"

namespace Collision
{

struct CollisionDetector
{
private:
    /** Responsible for building the BVH */
    LBVHBuilder _lbvh_builder;
    
    /** Responsible for traversing the BVH */
    LBVHTraversal _lbvh_traversal;

    /** Cache for potential collisions. */
    std::vector<std::pair<unsigned, unsigned>> _potential_collisions;

    /** Last iteration's detected collisions and accompanying order sorted by key */
    std::vector<DetectedCollision> _prev_detected_collisions;
    std::vector<unsigned> _prev_sorted_order;

    /** Current detected collisions and accompanying order sorted by key */
    std::vector<DetectedCollision> _cur_detected_collisions;
    std::vector<unsigned> _cur_sorted_order;

    /** Offsets for merging detected collisions after parallel narrow-phase */
    std::vector<unsigned> _merge_offsets;

    /** Generates a unique constexpr collision key for each pair of colliding objects that we can switch over. */
    constexpr static unsigned _makeCollisionKey(CollisionGeometryType a, CollisionGeometryType b)
    {
        unsigned ia = static_cast<unsigned>(a);
        unsigned ib = static_cast<unsigned>(b);

        if (ia > ib)
            std::swap(ia, ib);

        return (ia << 8) | ib;
    }
    
    /** Helper to skip adjacent primitives before running narrow-phase collision detection */
    inline bool _shouldSkip(const CollisionPrimitivePool& cpool, unsigned pi, unsigned pj);

    /** Perform narrow-phase collision detection. */
    inline void _narrowPhaseCollisionDetection(Sim::SimulationContext& ctx);
    
    /** Perform narrow-phase collision detection in parallel. */
    inline void _narrowPhaseCollisionDetection_Parallel(WorkerThreadContext& w_ctx, Sim::SimulationContext& ctx);

    /** Merge detected collisions from each thread into global detected collisions vector */
    inline void _mergeDetectedCollisions(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, std::vector<DetectedCollision>& merged_detected_collisions);

    /** Check a pair of primitives a and b for collision */
    inline void _processPotentialCollision(Sim::SimulationContext& ctx, unsigned a, unsigned b, std::vector<DetectedCollision>& detected_collisions);

    /** Specific subroutines for primitive-primitive narrow-phase collision checks */
    inline void _triangleSphere(Sim::SimulationContext& ctx, unsigned triangle, unsigned sphere, std::vector<DetectedCollision>& detected_collisions);
    inline void _triangleTriangle(Sim::SimulationContext& ctx, unsigned triangle1, unsigned triangle2, std::vector<DetectedCollision>& detected_collisions);
    inline void _rodRod(Sim::SimulationContext& ctx, unsigned rod1, unsigned rod2, std::vector<DetectedCollision>& detected_collisions);
    inline void _sphereSphere(Sim::SimulationContext& ctx, unsigned sphere1, unsigned sphere2, std::vector<DetectedCollision>& detected_collisions);

    /** General triangle-SDF continuous collision detection
     * Follows the implementation described by Pelletier-Guenette et al (2025): https://dl.acm.org/doi/full/10.1145/3747862
     * @param ctx simulation context
     * @param triangle primitive index for the triangle
     * @param rb primitive index for the rigid body
     * @param dt the time step size of the sim
     * @param normal (output) the collision normal
     * @param cp_barys (output) the barycentric coordinates of the contact point on the triangle
     * @param cp_rb_local (output) the contact point on the rigid body, expressed in the local rigid body frame
     * @returns whether or not a collision was detected over the interval [0, dt] given the current velocities
     */
    inline bool _triangleSDF_CCD(Sim::SimulationContext& ctx, unsigned triangle, unsigned rb, Real dt, Vec3r& normal, Vec3r& cp_barys, Vec3r& cp_rb_local);

    inline void _handleDetectedCollisions(Sim::SimulationContext& ctx);
    inline void _addCollision(Sim::SimulationContext& ctx, DetectedCollision& collision);
    inline void _removeCollision(Sim::SimulationContext& ctx, DetectedCollision& collision);
    inline void _updateCollision(Sim::SimulationContext& ctx, DetectedCollision& collision);
    

public:
    CollisionDetector() = default;

    /** Reserves memory for caches */
    CollisionDetector(unsigned capacity);

    /** Performs entire collision detection process: 
     * 
     * - Builds BVH 
     * - uses BVH for broad-phase collision detection
     * - Narrow-phase collision detection on potential collisions
     * - Creates appropriate collision energies/constraints and adds them to the energy pool.
     * - Recreates the adjacency graph and coloring.
     */
    void detectCollisionsAndRecolor(Sim::SimulationContext& ctx);

    /** Performs entire collision detection process in parallel */
    void detectCollisionsAndRecolor_Parallel(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, Sim::SimulationContext& ctx);
};

} // namespace Collision