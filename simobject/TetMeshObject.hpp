#pragma once

#include "simobject/ObjectBase.hpp"
#include "config/TetMeshObjectConfig.hpp"

#include "common/mesh/ParticleTetMesh.hpp"
#include "energy/TetElementEnergy.hpp"
#include "energy/GroundCollisionEnergy.hpp"

#include "energy/HardConstraintEnergy.hpp"
#include "energy/constraint/AttachmentConstraint.hpp"
#include "energy/constraint/GroundCollisionConstraint.hpp"

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

    /** Provides a way to iterate through all energies owned by the object. */
    virtual void for_each_energy(std::function<void(Energy::Energy_Base*)> func) override;
    virtual void for_each_energy(std::function<void(const Energy::Energy_Base*)> func) const override;

    /** Provides a way to iterate through all QUADRATIC energies owned by the object. */
    // virtual void for_each_quadratic_energy(std::function<void(QuadraticEnergy*)> func) override;
    // virtual void for_each_quadratic_energy(std::function<void(QuadraticEnergy*)> func) const override;

protected:
    std::string _filename;

    ParticleTetMesh _mesh;
    std::vector<Energy::TetElementEnergy> _element_energies;

    /** TODO: move this somewhere else. */
    std::vector<Energy::GroundCollisionEnergy> _collision_energies;

    /** TODO: temporary */
    std::vector<Energy::AttachmentConstraint> _attachment_constraints;
    std::vector<Energy::GroundCollisionConstraint> _ground_constraints;
    std::vector<Energy::HardConstraintEnergy> _constraint_energies;

    /** Material properties */
    Real _E;
    Real _nu;
    Real _density;

private:
    Config::TetMeshObjectConfig _config;
};

} // namespace SimObject