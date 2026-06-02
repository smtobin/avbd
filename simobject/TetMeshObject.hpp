#pragma once

#include "simobject/ObjectBase.hpp"
#include "config/TetMeshObjectConfig.hpp"

#include "common/mesh/ParticleTetMesh.hpp"
#include "energy/TetElementEnergy.hpp"
#include "energy/GroundCollisionEnergy.hpp"

namespace SimObject
{

class TetMeshObject : public Object_Base
{
public:
    TetMeshObject(const Config::TetMeshObjectConfig& config);
    
    const ParticleTetMesh* mesh() const { return &_mesh; }
    ParticleTetMesh* mesh() { return &_mesh; }

    virtual void setup() override;

    /** Provides a way to iterate through the particles owned by the object. */
    virtual void for_each_particle(std::function<void(Particle*)> func) override;
    virtual void for_each_particle(std::function<void(const Particle*)> func) const override;

protected:
    std::string _filename;

    ParticleTetMesh _mesh;
    std::vector<Energy::TetElementEnergy> _element_energies;

    /** TODO: move this somewhere else. */
    std::vector<Energy::GroundCollisionEnergy> _collision_energies;

private:
    Config::TetMeshObjectConfig _config;
};

} // namespace SimObject