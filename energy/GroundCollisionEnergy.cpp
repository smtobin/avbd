#include "energy/GroundCollisionEnergy.hpp"

namespace Energy
{

GroundCollisionEnergy::GroundCollisionEnergy(const Particle* particle)
    : QuadraticEnergy(1e10*ConstraintVecType::Ones()),
    _particle(particle)
{

}

const Particle* GroundCollisionEnergy::particle(int index) const
{
    return _particle;
}

GroundCollisionEnergy::ConstraintVecType GroundCollisionEnergy::evaluateConstraint() const
{
    ConstraintVecType C;
    C[0] = std::max(Real(0), -_particle->position[1]);
    return C;
}

/** Computes the gradient of the energy with respect to a particular particle. */
Vec3r GroundCollisionEnergy::gradient(int index, Real /* dt */) const
{
    if (_particle->position[1] >= 0)
    {
        _particle->in_collision = false;
        return Vec3r::Zero();
    }
    else
    {
        _particle->in_collision = true;
        return Vec3r(0, _k_cur[0]*_particle->position[1], 0);
    }
}

/** Computes the Hessian of the energy function with respect to a particular particle. */
Mat3r GroundCollisionEnergy::hessian(int index, Real /* dt */) const
{
    if (_particle->position[1] >= 0)
        return Mat3r::Zero();
    else
    {
        Mat3r hess = Mat3r::Zero();
        hess(1, 1) = _k_cur[0];
        return hess;
    }
}

} // namespace Energy