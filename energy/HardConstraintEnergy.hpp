#pragma once

#include "energy/EnergyBase.hpp"
#include "energy/constraint/ConstraintBase.hpp"

#include <numeric>

namespace Energy
{

class HardConstraintEnergy : public Energy_Base
{
public:
    HardConstraintEnergy(const Constraint_Base* constraint, Real k_start, 
        Real lambda_min=std::numeric_limits<Real>::lowest(), Real lambda_max=std::numeric_limits<Real>::max());

    /** Number of particles affected by the energy expression. */
    virtual int numParticles() const override { return _constraint->numParticles(); }

    /** The i'th particle affected by the energy expression. */
    virtual const Particle* particle(int index) const override { return _constraint->particle(index); };

    /** Resets the energy for a new time step */
    virtual void reset() override;

    /** Updates the energy after an iteration.
     * E.g. increases the stiffness, updates Lagrange multipliers, etc.
     */
    virtual void updateAfterIteration() override;

    /** Returns the current energy given the current state. */
    virtual Real energy(Real dt) const override;

    /** Computes the gradient of the energy with respect to a particular particle. */
    virtual Vec3r gradient(int index, Real dt) const override;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    virtual Mat3r hessian(int index, Real dt) const override; 

private:
    /** Pointer to the constraint
     * 
     * TODO: should this class own the constraint and be a templated class?
     */
    const Constraint_Base* _constraint;

    /** The current stiffness of the hard constraint. */
    Real _k;

    /** The initial stiffness. */
    Real _k_start;

    /** The current Lagrange multiplier. */
    Real _lambda;

    /** Bounds on the Lagrange multiplier (useful for inequality constraints) */
    Real _lambda_min;
    Real _lambda_max;
};

} // namespace Energy