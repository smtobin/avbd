#pragma once

#include "common/common.hpp"
#include "collision/CollisionPrimitivePool.hpp"
#include "collision/DetectedCollision.hpp"

namespace Collision
{

struct CollisionDetector
{
private:
    /** Cache for potential collisions. */
    std::vector<std::pair<unsigned, unsigned>> _potential_collisions;

    /** Last iteration's detected collisions and accompanying order sorted by key */
    std::vector<DetectedCollision> _prev_detected_collisions;
    std::vector<unsigned> _prev_sorted_order;

    /** Current detected collisions and accompanying order sorted by key */
    std::vector<DetectedCollision> _cur_detected_collisions;
    std::vector<unsigned> _cur_sorted_order;

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
    inline void _narrowPhaseCollisionDetection(Sim::SimulationContext& ctx, const std::vector<std::pair<unsigned, unsigned>>& potential_collisions);

    /** Specific subroutines for primitive-primitive narrow-phase collision checks */
    inline void _triangleSphere(Sim::SimulationContext& ctx, unsigned triangle, unsigned sphere);
    inline void _triangleTriangle(Sim::SimulationContext& ctx, unsigned triangle1, unsigned triangle2);
    inline void _sphereSphere(Sim::SimulationContext& ctx, unsigned sphere1, unsigned sphere2);
    

public:
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
};

} // namespace Collision