#pragma once

#include "config/Config.hpp"
#include "config/ObjectRenderConfig.hpp"


namespace Config
{

class ObjectConfig : public Config_Base
{
public:
    
    explicit ObjectConfig()
        : Config_Base(), _render_config()
    {}

    explicit ObjectConfig(const YAML::Node& node)
        : Config_Base(node), _render_config(node)
    {
        // _extractParameter("collisions", node, _collisions);
        // _extractParameter("graphics-only", node, _graphics_only);

        _extractParameter("position", node, _initial_position);
        _extractParameter("velocity", node, _initial_velocity);
        _extractParameter("rotation", node, _initial_rotation);
        _extractParameter("angular-velocity", node, _initial_angular_velocity);

        // _extractParameter("material", node, _material_class);
    }

    explicit ObjectConfig(
        const std::string& name,
        const Vec3r& initial_position,
        const Vec3r& initial_rotation,
        const Vec3r& initial_velocity,
        const Vec3r& initial_angular_velocity
    )
        : Config_Base(name), _render_config()
    {
        // _material_class.value = material_class;
        
        _initial_position.value = initial_position;
        _initial_rotation.value = initial_rotation;
        _initial_velocity.value = initial_velocity;
        _initial_angular_velocity.value = initial_angular_velocity;
        // _collisions.value = collisions;
        // _graphics_only.value = graphics_only;
    }
    
    // bool collisions() const { return _collisions.value; }
    // bool graphicsOnly() const { return _graphics_only.value; }
    Vec3r initialPosition() const { return _initial_position.value; }
    Vec3r initialVelocity() const { return _initial_velocity.value; }
    Vec3r initialRotation() const { return _initial_rotation.value; }
    Vec3r initialAngularVelocity() const { return _initial_angular_velocity.value; }

    // std::string materialClass() const { return _material_class.value; }

    const ObjectRenderConfig& renderConfig() const { return _render_config; }

protected:

    ConfigParameter<bool> _collisions = ConfigParameter<bool>(false);
    ConfigParameter<bool> _graphics_only = ConfigParameter<bool>(false);

    ConfigParameter<Vec3r> _initial_position = ConfigParameter<Vec3r>(Vec3r(0,0,0));
    ConfigParameter<Vec3r> _initial_velocity = ConfigParameter<Vec3r>(Vec3r(0,0,0));
    ConfigParameter<Vec3r> _initial_rotation = ConfigParameter<Vec3r>(Vec3r(0,0,0));
    ConfigParameter<Vec3r> _initial_angular_velocity = ConfigParameter<Vec3r>(Vec3r(0,0,0));

    ConfigParameter<std::string> _material_class = ConfigParameter<std::string>("default");

    ObjectRenderConfig _render_config;
};

} // namespace Config