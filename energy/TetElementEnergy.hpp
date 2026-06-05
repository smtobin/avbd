#pragma once

#include "energy/EnergyBase.hpp"
#include "common/Particle.hpp"
#include "common/mesh/ParticleTetMesh.hpp"

#include <array>

namespace Energy
{

class TetElementEnergy : public Energy_Base
{
public:
    TetElementEnergy(const ParticleTetMesh* mesh, int element_index, Real lambda, Real mu, Real kd);
    
    /** Number of particles affected by the energy expression. */
    virtual int numParticles() const override { return 4; }

    /** The i'th particle affected by the energy expression. */
    virtual const Particle* particle(int index) const override;

    /** Returns the current energy given the current state. */
    virtual Real energy(Real dt) const override;

    /** Computes the gradient of the energy with respect to a particular particle. */
    virtual Vec3r gradient(int index, Real dt) const override;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    virtual Mat3r hessian(int index, Real dt) const override; 

protected:
    /** Pointer to the mesh that owns this element */
    const ParticleTetMesh* _mesh;
    /** Element index in the tet mesh */
    int _element;

    /** Element material properties (Lame params) */
    Real _lambda;
    Real _mu;

    /** Strain-rate damping coefficient */
    Real _kd;
};


} // namespace Energy