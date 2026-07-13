#include "collision/CollisionObject.hpp"

#include "simobject/rigid/RigidSphere.hpp"
#include "simobject/TetMeshObject.hpp"

namespace Collision
{

CollisionObject::CollisionObject(SimObject::TetMeshObject* tet_mesh)
    : obj(tet_mesh)
    , type(ColliderType::TET_MESH)
    , fixed(false)
{}

CollisionObject::CollisionObject(SimObject::RigidSphere* sphere)
    : obj(sphere)
    , type(ColliderType::SPHERE)
    , fixed(sphere->fixed())
{}

} // namespace Collision