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
    return (_particle->position - _attachment_pos).norm();
}

Vec3r AttachmentConstraint::gradient(int /* index */) const
{
    Vec3r diff = _particle->position - _attachment_pos;
    Real norm = diff.norm();

    if (std::abs(norm) > 1e-12)
        return diff / norm;
    else
        return Vec3r(0,0,1);
}

Mat3r AttachmentConstraint::hessian(int /* index */) const
{
    Vec3r diff = _particle->position - _attachment_pos;
    Real norm = diff.norm();

    if (std::abs(norm) < 1e-12)
        return Mat3r::Identity();

    Mat3r hess = Mat3r::Identity() / norm - (diff * diff.transpose()) / (norm*norm*norm);
    return hess;
}
    
} // namespace Energy