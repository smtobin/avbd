#pragma once

#include "energy/constraint/ConstraintBase.hpp"
#include "common/Particle.hpp"

namespace Energy
{

/** Defines a constraint that fixes a particle to a particular position. */
class GroundCollisionConstraint : public Constraint_Base
{
public:
    GroundCollisionConstraint(const Particle* particle);
    virtual ~GroundCollisionConstraint() = default;

    /** Number of particles involved in this constraint. */
    virtual int numParticles() const { return 1; };

    /** The i'th particle affected by the energy expression. */
    virtual const Particle* particle(int /* index */) const override { return _particle; };

    /** Evaluates the constraint C(x) */
    virtual Real evaluate() const override;

    /** Evaluates the constraint gradient w.r.t. a particular particle */
    virtual Vec3r gradient(int index) const override;

    /** Evaluates the constraint Hessian w.r.t. a particular particle */
    virtual Mat3r hessian(int index) const override;


private:
    /** The particle being constrained */
    const Particle* _particle;
};

} // namespace Energy