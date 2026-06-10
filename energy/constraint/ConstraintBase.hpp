#pragma once

#include "common/common.hpp"
#include "common/Particle.hpp"

namespace Energy
{

/** Defines a common base class for constraints. */
class Constraint_Base
{
public:
    Constraint_Base() = default;
    virtual ~Constraint_Base() = default;

    /** Number of particles involved in this constraint. */
    virtual int numParticles() const = 0;

    /** The i'th particle affected by the constraint. */
    virtual const Particle* particle(int index) const = 0;

    /** Evaluates the constraint C(x) */
    virtual Real evaluate() const = 0;

    /** Evaluates the constraint gradient w.r.t. a particular particle */
    virtual Vec3r gradient(int index) const = 0;

    /** Evaluates the constraint Hessian w.r.t. a particular particle */
    virtual Mat3r hessian(int index) const = 0;
};

} // namespace Energy