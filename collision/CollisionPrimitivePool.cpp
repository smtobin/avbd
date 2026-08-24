#include "collision/CollisionPrimitivePool.hpp"

#include "simobject/TetMeshObject.hpp"
#include "simobject/Rod.hpp"
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

void CollisionPrimitivePool::addObject(const SimObject::Rod& rod)
{
    // create an entry in the SDFPrimitivePool to store the rod parameters required for collision
    unsigned sdf_slot = sdf_pool.addObject(rod);

    // add a rod segment primitive for each element in the rod
    // (assuming the rod is linear for now)
    const std::vector<unsigned>& nodes = rod.nodes();
    for (unsigned n_idx = 0; n_idx < nodes.size()-1; n_idx++)
    {
        unsigned slot = allocSlot();
        type[slot] = PrimitiveType::RodSegment;
        particle_indices[slot][0] = nodes[n_idx];
        particle_indices[slot][1] = nodes[n_idx+1];
        particle_indices[slot][2] = sdf_slot;   // index 2 = slot in the SDFPrimitivePool for auxiliary info
        num_particles[slot] = 2;
        object_id[slot] = rod.id();
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