#pragma once

#include "energy/collision/CollisionConstraintEnergyDefs.hpp"
#include "collision/CollisionShapeParams.hpp"

namespace Energy
{

/** Required data for triangle-rigid collision energy */
struct TriangleRigidCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec4u particle_indices;     // particle indices - first 3 indices are the triangle vertices, last index is the rigid body
    Collision::CollisionShapeParams* sdf_params;  // pointer to SDF parameters - will be used to reevaluate the constraint
    Vec3r barys;    // barycentric coordinates of the contact point on the face
    Vec3r cp_rb_local;  // contact point on the rigid body, in the rigid body's local frame
};

////////////////////////////////////////////////////////////////////////////////////////////

/** A non-owning view into the triangle-rigid collision energy storage */
struct TriangleRigidCollisionEnergyPoolView : CollisionConstraintEnergyPoolView<TriangleRigidCollisionEnergyInfo>
{
    using PoolType = TriangleRigidCollisionEnergyPool;
};

} // namespace Energy