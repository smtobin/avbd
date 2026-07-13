#pragma once

#include "config/RigidObjectConfig.hpp"

namespace SimObject
{
    class RigidSphere;
}

namespace Config
{

class RigidSphereConfig : public RigidObjectConfig
{
public:
    using SimObjectType = SimObject::RigidSphere;

    explicit RigidSphereConfig()
        : RigidObjectConfig()
    {}

    explicit RigidSphereConfig(const YAML::Node& node)
        : RigidObjectConfig(node)
    {
        _extractParameter("radius", node, _radius);
    }

    // explicit RigidSphereConfig(const std::string& name, const Vec3r& initial_position, const Vec3r& initial_rotation,
    //     const Vec3r& initial_velocity, const Vec3r& initial_angular_velocity, bool collisions,
    //     Real density, bool fixed,
    //     Real radius)
    //     : XPBDRigidBodyConfig(name, initial_position, initial_rotation, initial_velocity, initial_angular_velocity, collisions, density, fixed)
    // {
    //     _radius.value = radius;
    // }

    Real radius() const { return _radius.value; }

private:
    ConfigParameter<Real> _radius = ConfigParameter<Real>(0.5);

};

} // namespace Config