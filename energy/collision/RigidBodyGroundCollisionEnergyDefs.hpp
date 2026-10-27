#pragma once

#include "energy/collision/CollisionConstraintEnergyDefs.hpp"
#include "collision/CollisionShapeParams.hpp"

namespace Energy
{

/** Required data for rigid body ground collision energies. */
struct RigidBodyGroundCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec1u particle_indices;     // particle indices for each constraint
    Vec3r cp_rb_local;  // local vector to contact point on rigid body
    Vec3r cp_ground;    // contact point on the ground
    Collision::CollisionShapeParams shape_params;     // describes the rigid body shape
};

/////////////////////////////////////////////////////////////////////////////////////////

/** A non-owning view into the rigid body ground collision energy storage. */
struct RigidBodyGroundCollisionEnergyPoolView : CollisionConstraintEnergyPoolView<RigidBodyGroundCollisionEnergyInfo>
{
};

} // namespace Energy