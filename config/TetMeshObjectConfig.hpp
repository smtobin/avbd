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
    }

    std::string filename() const { return _filename.value; }

    const MeshRenderConfig& meshRenderConfig() const { return _mesh_render_config; }

protected:
    ConfigParameter<std::string> _filename = ConfigParameter<std::string>("");

    MeshRenderConfig _mesh_render_config;
};

} // namespace Config