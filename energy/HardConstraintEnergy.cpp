#include "energy/HardConstraintEnergy.hpp"

namespace Energy
{

HardConstraintEnergy::HardConstraintEnergy(const Constraint_Base* constraint, Real k_start)
    : Energy_Base(),
    _constraint(constraint),
    _k(0), _k_start(k_start),
    _lambda(0)
{

}

void HardConstraintEnergy::reset()
{
    _k = std::max(STIFFNESS_GAMMA*_k, _k_start);
    _lambda = STIFFNESS_GAMMA * _lambda;

    /** TODO: incorporate alpha */
}

void HardConstraintEnergy::updateAfterIteration()
{
    Real C = _constraint->evaluate();

    _lambda += _k * C;
    _k += STIFFNESS_BETA * C;
}

Real HardConstraintEnergy::energy(Real /* dt */) const
{
    Real C = _constraint->evaluate();
    return 0.5 * _k * C * C + _lambda * C;
}

Vec3r HardConstraintEnergy::gradient(int index, Real /* dt */) const
{
    Real C = _constraint->evaluate();
    Vec3r gradC = _constraint->gradient(index);

    return (_k * C + _lambda) * gradC;
}

Mat3r HardConstraintEnergy::hessian(int index, Real /* dt */) const
{
    Real C = _constraint->evaluate();
    Vec3r gradC = _constraint->gradient(index);
    Mat3r hessC = _constraint->hessian(index);

    return (_k * C + _lambda) * hessC + _k * gradC * gradC.transpose();
}

} // namespace Energy