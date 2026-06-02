#include "energy/GroundCollisionEnergy.hpp"

namespace Energy
{

GroundCollisionEnergy::GroundCollisionEnergy(const Particle* particle)
    : Energy_Base(),
    _particle(particle), _k(1e10)
{

}

const Particle* GroundCollisionEnergy::particle(int index) const
{
    return _particle;
}

/** Returns the current energy given the current state. */
Real GroundCollisionEnergy::energy() const
{
    Real dist = std::max(Real(0), -_particle->position[1]);
    return 0.5*_k*dist*dist;
}

/** Computes the gradient of the energy with respect to a particular particle. */
Vec3r GroundCollisionEnergy::gradient(int index) const
{
    if (_particle->position[1] >= 0)
        return Vec3r::Zero();
    else
    {
        return Vec3r(0, -0.5*_k*_particle->position[1], 0);
    }
}

/** Computes the Hessian of the energy function with respect to a particular particle. */
Mat3r GroundCollisionEnergy::hessian(int index) const
{
    if (_particle->position[1] >= 0)
        return Mat3r::Zero();
    else
    {
        Mat3r hess = Mat3r::Zero();
        hess(1, 1) = -0.5*_k;
        return hess;
    }
}

} // namespace Energy