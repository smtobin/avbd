#pragma once

#include "collision/CollisionGeometryType.hpp"

namespace Collision
{

/** Shape parameters required for collision detection, e.g. radius, height, etc.
 * Used primarily for analytic rigid shapes.
*/
struct CollisionShapeParams
{
    CollisionGeometryType type;
    union {
        struct { Real radius; } sphere;
        struct { Real half_extents[3]; } box;
        struct { Real radius, half_height; } capsule;
        struct { Real radius; } rod;
    };
};

} // namespace Collision