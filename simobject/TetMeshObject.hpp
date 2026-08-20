#pragma once

#include "simobject/ObjectBase.hpp"
#include "config/TetMeshObjectConfig.hpp"

#include "common/mesh/ParticleTetMesh.hpp"

namespace SimObject
{

class TetMeshObject : public Object_Base
{
public:
    TetMeshObject(Sim::SimulationContext* ctx, const Config::TetMeshObjectConfig& config);
    
    const ParticleTetMesh& mesh() const { return _mesh; }
    ParticleTetMesh& mesh() { return _mesh; }

    virtual void setup() override;

protected:
    std::string _filename;

    ParticleTetMesh _mesh;

    /** Material properties */
    Real _E;
    Real _nu;
    Real _density;

private:
    Config::TetMeshObjectConfig _config;
};

} // namespace SimObject