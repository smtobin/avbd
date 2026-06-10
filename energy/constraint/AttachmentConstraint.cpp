#include "energy/constraint/AttachmentConstraint.hpp"

namespace Energy
{

AttachmentConstraint::AttachmentConstraint(const Particle* particle, const Vec3r& attachment_pos)
    : Constraint_Base(),
     _particle(particle), _attachment_pos(attachment_pos)
{

}

Real AttachmentConstraint::evaluate() const
{
    return 0.5*(_particle->position - _attachment_pos).squaredNorm();
}

Vec3r AttachmentConstraint::gradient(int /* index */) const
{
    return _particle->position - _attachment_pos;
}

Mat3r AttachmentConstraint::hessian(int /* index */) const
{
    return Mat3r::Identity();
}
    
} // namespace Energy