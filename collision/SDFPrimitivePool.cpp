#include "collision/SDFPrimitivePool.hpp"

#include "simobject/rigid/RigidSphere.hpp"

namespace Collision
{

unsigned SDFPrimitivePool::addObject(const SimObject::RigidSphere& sphere)
{
    unsigned slot = allocSlot();
    params[slot] = {
        SDFType::Sphere,
        { sphere.radius() }
    };
    particles[slot] = sphere.com();

    return slot;
}

} // namespace Collision