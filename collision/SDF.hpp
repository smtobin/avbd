#pragma once

#include "common/common.hpp"
#include "collision/AABB.hpp"

namespace Collision
{

/** Enum for the different supported SDFs
 * 
 * Sphere : sphere centered at the origin
 * Box : box centered at the origin
 * Capsule : capsule centered at the origin, "height" direction along y-axis
 */
enum class SDFType : uint8_t {
    Sphere,
    Box,
    Capsule
};

struct SDFShapeParams
{
    SDFType type;
    union {
        struct { Real radius; } sphere;
        struct { Real half_extents[3]; } box;
        struct { Real radius, half_height; } capsule;
    };
};

struct SDF
{

    /** Returns the AABB for the object in local space based on its type
     * @param s : SDF parameters
     */
    static inline AABB localBounds(const SDFShapeParams& s)
    {
        switch (s.type)
        {
            case SDFType::Sphere:
            {
                return { Vec3r::Constant(-s.sphere.radius), Vec3r::Constant(s.sphere.radius) };
            }
            case SDFType::Box:
            {
                return { -Vec3r(s.box.half_extents), Vec3r(s.box.half_extents) };
            }
            case SDFType::Capsule:
            {  
                return { 
                    Vec3r(-s.capsule.radius, -s.capsule.half_height - s.capsule.radius, -s.capsule.radius),
                    Vec3r( s.capsule.radius,  s.capsule.half_height + s.capsule.radius, s.capsule.radius) 
                };
            }
        }
    }

    /** Evaluates the SDF based on its type
     * @param s : SDF parameters
     * @param p : the query point (in local frame)
     */
    static inline Real evaluate(const SDFShapeParams& s, const Vec3r& p)
    {
        switch (s.type)
        {
            case SDFType::Sphere:
            {
                return p.norm() - s.sphere.radius;
            }
            /** TODO: (07/19/26) implement box SDF */
            case SDFType::Box:
            {
                throw std::runtime_error("Box SDF evaluate  not implemented yet.");
                return 0.0;
            }
            /** TODO: (07/19/26) implement capsule SDF */
            case SDFType::Capsule:
            {
                throw std::runtime_error("Capsule SDF evaluate not implemented yet.");  
                return 0.0;
            }
        }
    }

    /** Evaluates the SDF gradient based on its type
     * @param s : SDF parameters
     * @param p : the query point (in local frame)
     */
    static inline Vec3r gradient(const SDFShapeParams& s, const Vec3r& p)
    {
        switch (s.type)
        {
            case SDFType::Sphere:
            {
                Real dist = p.norm();
                if (dist > 1e-8)
                    return p / dist;
                else
                    return Vec3r(1,0,0);
            }
                
            /** TODO: (07/19/26) implement box SDF gradient */
            case SDFType::Box:
            {
                throw std::runtime_error("Box SDF gradient not implemented yet.");
            }
            /** TODO: (07/19/26) implement capsule SDF */
            case SDFType::Capsule:
            {
                throw std::runtime_error("Capsule SDF gradient not implemented yet.");  
            }
        }
    }

};

} // namespace Collision