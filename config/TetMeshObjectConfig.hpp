#pragma once

#include "config/ObjectConfig.hpp"
#include "config/MeshRenderConfig.hpp"

namespace Config
{

class TetMeshObjectConfig : public ObjectConfig
{
public:
    using SimObjectType = SimObject::TetMeshObject;

    explicit TetMeshObjectConfig()
        : ObjectConfig()
    {}

    explicit TetMeshObjectConfig(const YAML::Node& node)
        : ObjectConfig(node), _mesh_render_config(node)
    {
        _extractParameter("filename", node, _filename);

        _extractParameter("E", node, _E);
        _extractParameter("nu", node, _nu);
        _extractParameter("density", node, _density);

        _extractParameter("damping-coeff", node, _damping_coeff);

        _extractParameter("scaling", node, _scaling);
    }

    std::string filename() const { return _filename.value; }

    Vec3r scaling() const { return _scaling.value; }
    Real E() const { return _E.value; }
    Real nu() const { return _nu.value; }
    Real density() const { return _density.value; }

    Real dampingCoefficient() const { return _damping_coeff.value; }

    const MeshRenderConfig& meshRenderConfig() const { return _mesh_render_config; }

protected:
    ConfigParameter<std::string> _filename = ConfigParameter<std::string>("");

    ConfigParameter<Vec3r> _scaling = ConfigParameter<Vec3r>(Vec3r(1,1,1));

    ConfigParameter<Real> _E = ConfigParameter<Real>(1e6);
    ConfigParameter<Real> _nu = ConfigParameter<Real>(0.4);
    ConfigParameter<Real> _density = ConfigParameter<Real>(1000);

    ConfigParameter<Real> _damping_coeff = ConfigParameter<Real>(0);

    MeshRenderConfig _mesh_render_config;
};

} // namespace Config