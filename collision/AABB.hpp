#pragma once

#include "common/common.hpp"

namespace Collision
{

/** Axis-Aligned Bounding Box (AABB) */
struct AABB
{
    Vec3r min;
    Vec3r max;

    static AABB empty()
    {
        return { Vec3r(std::numeric_limits<Real>::max()), Vec3r(std::numeric_limits<Real>::lowest()) };
    }

    void expand(const Vec3r& p)
    {
        min = min.cwiseMin(p);
        max = max.cwiseMax(p);
    }

    void expand(const AABB& b)
    {
        min = min.cwiseMin(b.min);
        max = max.cwiseMax(b.max);
    }

    void pad(Real margin)
    {
        min -= margin;
        max += margin;
    }

    Vec3r center()  const { return (min + max) * Real(0.5); }
    Vec3r extent()  const { return (max - min) * Real(0.5); }

    bool overlaps(const AABB& o) const
    {
        // return (min.array() <= o.max.array()).all() &&
        //        (max.array() >= o.min.array()).all();

        return  (min[0] <= o.max[0]) && (min[1] <= o.max[1]) && (min[2] <= o.max[2]) &&
                (max[0] >= o.min[0]) && (max[1] >= o.min[1]) && (max[2] >= o.min[2]);
    }
};


inline std::ostream& operator<<(std::ostream& os, const Collision::AABB& aabb)
{
    os << aabb.min << " to " << aabb.max;
    return os;
}

} // namespace Collision