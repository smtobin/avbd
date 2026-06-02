#pragma once

#include "energy/QuadraticEnergy.hpp"
#include "common/Particle.hpp"

#include <array>

namespace Energy
{

/** Energy associated with particle collision with ground.
 * Quadratic energy approximation, assumes ground is y=0
 */
class GroundCollisionEnergy : public QuadraticEnergy<1>
{
public:
    using QuadraticEnergy<1>::ConstraintVecType;

    GroundCollisionEnergy(const Particle* particle);
    
    /** Number of particles affected by the energy expression. */
    virtual int numParticles() const override { return 1; }

    /** The i'th particle affected by the energy expression. */
    virtual const Particle* particle(int index) const override;

    /** Computes the gradient of the energy with respect to a particular particle. */
    virtual Vec3r gradient(int index) const override;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    virtual Mat3r hessian(int index) const override; 

    /** Evaluate C(x) */
    virtual ConstraintVecType evaluateConstraint() const override;

protected:
    /** Pointer to the particle */
    const Particle* _particle;
};


} // namespace Energy