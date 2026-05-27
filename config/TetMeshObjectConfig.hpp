#pragma once

#include "config/ObjectConfig.hpp"

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
        : ObjectConfig(node)
    {
        _extractParameter("filename", node, _filename);
    }

    std::string filename() const { return _filename.value; }

protected:
    ConfigParameter<std::string> _filename = ConfigParameter<std::string>("");
};

} // namespace Config