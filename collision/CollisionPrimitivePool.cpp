#include "collision/CollisionPrimitivePool.hpp"

#include "simobject/TetMeshObject.hpp"
#include "simobject/rigid/RigidSphere.hpp"

namespace Collision
{

void CollisionPrimitivePool::addObject(const SimObject::TetMeshObject& mesh_obj)
{
    const ParticleTetMesh& mesh = mesh_obj.mesh();
    const auto& vertices = mesh.vertices();
    for (const auto& f : mesh.faces())
    {
        unsigned slot = allocSlot();
        type[slot] = PrimitiveType::Triangle;
        particle_indices[slot][0] = vertices.at(f[0]);
        particle_indices[slot][1] = vertices.at(f[1]);
        particle_indices[slot][2] = vertices.at(f[2]);
        num_particles[slot] = 3;
        object_id[slot] = mesh_obj.id();
    }
}

unsigned CollisionPrimitivePool::addObject(const SimObject::RigidSphere& sphere)
{
    unsigned slot = allocSlot();
    type[slot] = PrimitiveType::RigidSDF;
    unsigned sdf_slot = sdf_pool.addObject(sphere);
    particle_indices[slot] = {sdf_slot};
    num_particles[slot] = 1;
    object_id[slot] = sphere.id();

    return slot;
}

} // namespace Collision