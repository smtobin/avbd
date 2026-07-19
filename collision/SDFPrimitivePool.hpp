#pragma once

#include "common/common.hpp"
#include "common/TombstonePool.hpp"
#include "collision/SDF.hpp"

namespace Collision
{

struct SDFPrimitivePool : TombstonePool
{
    std::vector<SDFShapeParams> params;     // SDF parameters
    std::vector<unsigned> particles;        // oriented particle indices - one for each center of mass

    /** Constructor initializes memory
     * @param capacity : the capacity of the memory pool
     */
    explicit SDFPrimitivePool(unsigned capacity)
        : TombstonePool(capacity)
        , params(capacity)
        , particles(capacity)
    {

    }

    unsigned addObject(const SimObject::RigidSphere& sphere);
};

} // namespace Collision