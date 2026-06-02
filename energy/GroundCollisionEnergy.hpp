#pragma once

#include "energy/EnergyBase.hpp"
#include "common/Particle.hpp"

#include <array>

namespace Energy
{

/** Energy associated with particle collision with ground.
 * Quadratic energy approximation, assumes ground is y=0
 */
class GroundCollisionEnergy : public Energy_Base
{
public:
    GroundCollisionEnergy(const Particle* particle);
    
    /** Number of particles affected by the energy expression. */
    virtual int numParticles() const override { return 1; }

    /** The i'th particle affected by the energy expression. */
    virtual const Particle* particle(int index) const override;

    /** Returns the current energy given the current state. */
    virtual Real energy() const override;

    /** Computes the gradient of the energy with respect to a particular particle. */
    virtual Vec3r gradient(int index) const override;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    virtual Mat3r hessian(int index) const override; 

protected:
    /** Pointer to the particle */
    const Particle* _particle;

    /** Energy stiffness */
    Real _k;
};


} // namespace Energy