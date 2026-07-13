#pragma once

#include "common/common.hpp"

namespace Collision
{

enum class ColliderType : uint8_t
{
    TET_MESH, SPHERE,
    COUNT
};

/** An object being tested for collision.
 * The object is stored as a type-tagged pointer, with type according to the ColliderType enum.
 */
struct CollisionObject
{
    void* obj;
    ColliderType type;
    bool fixed;

    CollisionObject(SimObject::TetMeshObject* tet_mesh);
    CollisionObject(SimObject::RigidSphere* sphere);
};



} // namespace Collision