#include "energy/constraint/GroundCollisionConstraint.hpp"

namespace Energy
{

GroundCollisionConstraint::GroundCollisionConstraint(const Particle* particle)
    : Constraint_Base(),
    _particle(particle)
{

}

Real GroundCollisionConstraint::evaluate() const
{
    if (_particle->position[1] <= 0)
    {
        _particle->in_collision = true;
    }
    
    return -_particle->position[1];
}

Vec3r GroundCollisionConstraint::gradient(int /* index */) const
{
    return Vec3r(0,-1,0);
}

Mat3r GroundCollisionConstraint::hessian(int /* index */) const
{
    return Mat3r::Zero();
}

} // namespace Energy