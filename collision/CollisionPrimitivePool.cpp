#include "collision/CollisionPrimitivePool.hpp"

#include "simobject/TetMeshObject.hpp"
#include "simobject/rigid/RigidSphere.hpp"

namespace Collision
{

void CollisionPrimitivePool::addObject(const SimObject::TetMeshObject& mesh_obj)
{
    const ParticleTetMesh& mesh = mesh_obj.mesh();
    for (const auto& f : mesh.faces())
    {
        unsigned slot = allocSlot();
        type[slot] = PrimitiveType::Triangle;
        particle_indices[slot][0] = f[0];
        particle_indices[slot][1] = f[1];
        particle_indices[slot][2] = f[2];
        num_particles[slot] = 3;
        object_id[slot] = mesh_obj.id();
    }
}

void CollisionPrimitivePool::addObject(const SimObject::RigidSphere& sphere)
{
    unsigned slot = allocSlot();
    type[slot] = PrimitiveType::RigidSDF;
    unsigned sdf_slot = sdf_pool.addObject(sphere);
    particle_indices[slot] = {sdf_slot};
    num_particles[slot] = 1;
    object_id[slot] = sphere.id();
}

} // namespace Collision