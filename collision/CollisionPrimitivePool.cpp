#include "collision/CollisionPrimitivePool.hpp"

#include "simobject/TetMeshObject.hpp"
#include "simobject/Rod.hpp"
#include "simobject/rigid/RigidSphere.hpp"

namespace Collision
{

unsigned CollisionShapeParamsPool::addObject(const SimObject::RigidSphere& sphere)
{
    unsigned slot = allocSlot();
    shape_params[slot] = {
        CollisionGeometryType::Sphere,
        { sphere.radius() }
    };

    return slot;
}

unsigned CollisionShapeParamsPool::addObject(const SimObject::Rod& rod)
{
    unsigned slot = allocSlot();
    shape_params[slot] = {
        CollisionGeometryType::RodSegment,
        { rod.radius() }
    };

    return slot;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


void CollisionPrimitivePool::addObject(const SimObject::TetMeshObject& mesh_obj)
{
    const ParticleTetMesh& mesh = mesh_obj.mesh();
    const auto& vertices = mesh.vertices();
    for (const auto& f : mesh.faces())
    {
        unsigned slot = allocSlot();
        type[slot] = CollisionGeometryType::Triangle;
        particle_indices[slot][0] = vertices.at(f[0]);
        particle_indices[slot][1] = vertices.at(f[1]);
        particle_indices[slot][2] = vertices.at(f[2]);
        num_particles[slot] = 3;
        object_id[slot] = mesh_obj.id();
    }
}

unsigned CollisionPrimitivePool::addObject(const SimObject::Rod& rod)
{
    // create an entry in the CollisionShapeParamsPool to store the rod parameters required for collision
    unsigned shape_params_slot = shape_params_pool.addObject(rod);

    // add a rod segment primitive for each element in the rod
    // (assuming the rod is linear for now)
    const std::vector<unsigned>& nodes = rod.nodes();
    for (unsigned n_idx = 0; n_idx < nodes.size()-1; n_idx++)
    {
        unsigned slot = allocSlot();
        type[slot] = CollisionGeometryType::RodSegment;
        particle_indices[slot][0] = nodes[n_idx];
        particle_indices[slot][1] = nodes[n_idx+1];
        particle_indices[slot][2] = shape_params_slot;   // index 2 = slot in the CollisionShapeParamsPool for auxiliary info
        num_particles[slot] = 2;
        object_id[slot] = rod.id();
    }

    return shape_params_slot;
}

unsigned CollisionPrimitivePool::addObject(const SimObject::RigidSphere& sphere)
{
    // create an entry in the CollisionShapeParamsPool to store sphere parameters required for collision
    unsigned shape_params_slot = shape_params_pool.addObject(sphere);

    unsigned slot = allocSlot();
    type[slot] = CollisionGeometryType::Sphere;
    particle_indices[slot][0] = sphere.com();
    particle_indices[slot][1] = shape_params_slot;
    num_particles[slot] = 1;
    object_id[slot] = sphere.id();

    return shape_params_slot;
}

} // namespace Collision