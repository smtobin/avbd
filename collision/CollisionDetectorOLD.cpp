#include "collision/CollisionDetector.hpp"

namespace Collision
{

// initialize static values
bool CollisionDetector::_collision_table_initialized = false;
CollisionDetector::CollisionFunc CollisionDetector::_collision_table[static_cast<int>(ColliderType::COUNT)][static_cast<int>(ColliderType::COUNT)] = {};

// [](void* a, void* b) {
    // CollisionScene::_checkCollision(static_cast<SimObject::XPBDRigidSphere*>(a), static_cast<SimObject::XPBDRigidSphere*>(b));
// };
void CollisionDetector::_initCollisionTable()
{
    if (_collision_table_initialized)
        return;
    
    // first type is a tet mesh
    _collision_table[static_cast<int>(ColliderType::TET_MESH)][static_cast<int>(ColliderType::TET_MESH)] = [](CollisionDetector* detector, void* a, void* b) {
        CollisionDetector::_checkCollision(detector, static_cast<SimObject::TetMeshObject*>(a), static_cast<SimObject::TetMeshObject*>(b));
    };
    _collision_table[static_cast<int>(ColliderType::TET_MESH)][static_cast<int>(ColliderType::SPHERE)] = [](CollisionDetector* detector, void* a, void* b) {
        CollisionDetector::_checkCollision(detector, static_cast<SimObject::TetMeshObject*>(a), static_cast<SimObject::TetMeshObject*>(b));
    };

    // first type is a sphere
    _collision_table[static_cast<int>(ColliderType::SPHERE)][static_cast<int>(ColliderType::SPHERE)] = [](CollisionDetector* detector, void* a, void* b) {
        CollisionDetector::_checkCollision(detector, static_cast<SimObject::RigidSphere*>(a), static_cast<SimObject::RigidSphere*>(b));
    };

    _collision_table[static_cast<int>(ColliderType::SPHERE)][static_cast<int>(ColliderType::TET_MESH)] = [](CollisionDetector* detector, void* a, void* b) {
        CollisionDetector::_checkCollision(detector, static_cast<SimObject::TetMeshObject*>(b), static_cast<SimObject::RigidSphere*>(a));
    };

    

    _collision_table_initialized = true;
}

CollisionDetector::CollisionDetector()
{

}

void CollisionDetector::_checkCollision(CollisionDetector* detector, SimObject::RigidSphere* sphere1, SimObject::RigidSphere* sphere2)
{
    // no sphere-sphere contact
    return;
}
void CollisionDetector::_checkCollision(CollisionDetector* detector, SimObject::TetMeshObject* mesh_obj, SimObject::RigidSphere* sphere)
{
    std::cout << "Checking for potential mesh-sphere collision" << std::endl;
}
void CollisionDetector::_checkCollision(CollisionDetector* detector, SimObject::TetMeshObject* mesh_obj1, SimObject::TetMeshObject* mesh_obj2)
{
    // no mesh-mesh contact
    return;
}



} // namespace Collision