#pragma once

#include "common/common.hpp"
#include "collision/CollisionShapeParams.hpp"
#include "collision/AABB.hpp"

namespace Collision
{

struct SDF
{

    /** Returns the AABB for the object in local space based on its type
     * @param s : SDF parameters
     */
    static inline AABB localBounds(const CollisionShapeParams& s)
    {
        switch (s.type)
        {
            case CollisionGeometryType::Sphere:
            {
                return { Vec3r::Constant(-s.sphere.radius), Vec3r::Constant(s.sphere.radius) };
            }
            case CollisionGeometryType::Box:
            {
                return { -Vec3r(s.box.half_extents), Vec3r(s.box.half_extents) };
            }
            case CollisionGeometryType::Capsule:
            {  
                return { 
                    Vec3r(-s.capsule.radius, -s.capsule.half_height - s.capsule.radius, -s.capsule.radius),
                    Vec3r( s.capsule.radius,  s.capsule.half_height + s.capsule.radius, s.capsule.radius) 
                };
            }
            default:
            {
                throw std::runtime_error("SDF::localBounds only implemented for analytic shape types.");
            }
        }
    }

    /** Evaluates the SDF based on its type
     * @param s : SDF parameters
     * @param p : the query point (in local frame)
     */
    static inline Real evaluate(const CollisionShapeParams& s, const Vec3r& p)
    {
        switch (s.type)
        {
            case CollisionGeometryType::Sphere:
            {
                return p.norm() - s.sphere.radius;
            }
            /** TODO: (07/19/26) implement box SDF */
            case CollisionGeometryType::Box:
            {
                throw std::runtime_error("Box SDF evaluate  not implemented yet.");
                return 0.0;
            }
            /** TODO: (07/19/26) implement capsule SDF */
            case CollisionGeometryType::Capsule:
            {
                throw std::runtime_error("Capsule SDF evaluate not implemented yet.");  
                return 0.0;
            }

            default:
            {
                throw std::runtime_error("SDF::evaluate oly implemented for analytic shape types.");
            }
        }
    }

    /** Evaluates the SDF gradient based on its type
     * @param s : SDF parameters
     * @param p : the query point (in local frame)
     */
    static inline Vec3r gradient(const CollisionShapeParams& s, const Vec3r& p)
    {
        switch (s.type)
        {
            case CollisionGeometryType::Sphere:
            {
                Real dist = p.norm();
                if (dist > 1e-8)
                    return p / dist;
                else
                    return Vec3r(1,0,0);
            }
                
            /** TODO: (07/19/26) implement box SDF gradient */
            case CollisionGeometryType::Box:
            {
                throw std::runtime_error("Box SDF gradient not implemented yet.");
            }
            /** TODO: (07/19/26) implement capsule SDF */
            case CollisionGeometryType::Capsule:
            {
                throw std::runtime_error("Capsule SDF gradient not implemented yet.");  
            }

            default:
            {
                throw std::runtime_error("SDF::gradient only implemented for analytic shape types.");
            }
        }
    }

};

} // namespace Collision