#pragma once

#include "common/common.hpp"
#include "collision/CollisionObject.hpp"

#include <vector>

namespace Collision
{

class CollisionDetector
{
private:
    /** Type-tagged generic collision objects */
    std::vector<CollisionObject> _collision_objects;

    /** Static dispatch table for colliding pairs of objects during narrow-phase collision detection */
    using CollisionFunc = void(*)(CollisionDetector* detector, void*, void*);
    static CollisionFunc _collision_table[static_cast<int>(ColliderType::COUNT)][static_cast<int>(ColliderType::COUNT)];

    /** Whether or not the collision table has been initialized yet */
    static bool _collision_table_initialized;

public:
    CollisionDetector();

    template <typename ObjectType>
    void addObject(ObjectType* obj)
    {
        auto& new_collision_obj = _collision_objects.emplace_back(obj);
    }

private:
    void _initCollisionTable();

    static void _checkCollision(CollisionDetector* detector, SimObject::RigidSphere* sphere1, SimObject::RigidSphere* sphere2);
    static void _checkCollision(CollisionDetector* detector, SimObject::TetMeshObject* mesh_obj, SimObject::RigidSphere* sphere);
    static void _checkCollision(CollisionDetector* detector, SimObject::TetMeshObject* mesh_obj1, SimObject::TetMeshObject* mesh_obj2);
};

} // Collision