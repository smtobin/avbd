#pragma once

#include "common/common.hpp"

namespace Collision
{

struct CollisionDetector
{
private:
    static bool _collision_table_initialized;
    // using NarrowPhaseFn = bool(*)(const Primitive&, const Primitive&, Contact&);
    // constexpr static NarrowPhaseFn _collision_table[(size_t)PrimitiveType::Count][(size_t)PrimitiveType::Count] =
    // {

    // }
    static void _initCollisionTable();

    

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