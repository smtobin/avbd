#pragma once

#include "common/common.hpp"

namespace Collision
{

struct CollisionDetector
{
private:
    /** Generates a unique constexpr collision key for each pair of colliding objects that we can switch over. */
    static constexpr unsigned _makeCollisionKey(CollisionGeometryType a, CollisionGeometryType b)
    {
        unsigned ia = static_cast<unsigned>(a);
        unsigned ib = static_cast<unsigned>(b);

        if (ia > ib)
            std::swap(ia, ib);

        return (ia << 8) | ib;
    }
    
    /** Helper to skip adjacent primitives before running narrow-phase collision detection */
    inline static bool _shouldSkip(const CollisionPrimitivePool& cpool, unsigned pi, unsigned pj);

    /** Perform narrow-phase collision detection. */
    inline static void _narrowPhaseCollisionDetection(const std::vector<std::pair<unsigned, unsigned>>& potential_collisions);

    /** Specific subroutines for primitive-primitive narrow-phase collision checks */
    inline static void _triangleSphere(Sim::SimulationContext& ctx, unsigned triangle, unsigned sphere);
    inline static void _triangleTriangle(Sim::SimulationContext& ctx, unsigned triangle1, unsigned triangle2);
    inline static void _sphereSphere(Sim::SimulationContext& ctx, unsigned sphere1, unsigned sphere2);
    

public:
    /** Performs entire collision detection process: 
     * 
     * - Builds BVH 
     * - uses BVH for broad-phase collision detection
     * - Narrow-phase collision detection on potential collisions
     * - Creates appropriate collision energies/constraints and adds them to the energy pool.
     * - Recreates the adjacency graph and coloring.
     */
    static void detectCollisionsAndRecolor(Sim::SimulationContext& ctx);
};

} // namespace Collision