#include "energy/HardConstraintEnergy.hpp"

namespace Energy
{

HardConstraintEnergy::HardConstraintEnergy(const Constraint_Base* constraint, Real k_start, Real lambda_min, Real lambda_max)
    : Energy_Base(),
    _constraint(constraint), _C_prev(0),
    _k(k_start), _k_start(k_start),
    _lambda(0), _lambda_min(lambda_min), _lambda_max(lambda_max)
{

}

Real HardConstraintEnergy::_evaluateConstraint() const
{
    return _constraint->evaluate() - CONSTRAINT_ALPHA*_C_prev;
}

void HardConstraintEnergy::reset()
{
    _k = std::max(STIFFNESS_GAMMA*_k, _k_start);
    _lambda = STIFFNESS_GAMMA * CONSTRAINT_ALPHA * _lambda;

    /** TODO: incorporate alpha */
}

void HardConstraintEnergy::updateAfterIteration()
{
    Real C = _evaluateConstraint();

    // std::cout << "HardConstraintEnergy C: " << C << std::endl;
    // std::cout << "  Lambda: " << _lambda << std::endl;
    // std::cout << "  k: " << _k << std::endl;

    Real lambda_p = _k * C + _lambda;
    _lambda = std::max(_lambda_min, std::min(_lambda_max, lambda_p));

    if (lambda_p > _lambda_min && lambda_p < _lambda_max)
    {
        _k += STIFFNESS_BETA * C;
        _C_prev = _constraint->evaluate();
    }

    
}

Real HardConstraintEnergy::energy(Real /* dt */) const
{
    Real C = _evaluateConstraint();
    return 0.5 * _k * C * C + _lambda * C;
}

Vec3r HardConstraintEnergy::gradient(int index, Real /* dt */) const
{
    Real C = _evaluateConstraint();
    Vec3r gradC = _constraint->gradient(index);

    Real lambda_p = std::max(_lambda_min, std::min(_lambda_max, _k*C + _lambda));

    return lambda_p * gradC;
}

Mat3r HardConstraintEnergy::hessian(int index, Real /* dt */) const
{
    Real C = _evaluateConstraint();
    Vec3r gradC = _constraint->gradient(index);
    Mat3r hessC = _constraint->hessian(index);

    Real lambda_p = _k*C + _lambda;

    Real k_scaled = _k;
    if (lambda_p < _lambda_min && std::abs(C) > 1e-12)
        k_scaled = (_lambda_min - _lambda) / C;
    else if (lambda_p > _lambda_max && std::abs(C) > 1e-12)
        k_scaled = (_lambda_max - _lambda) / C;

    Mat3r hess = (k_scaled * C + _lambda) * hessC + k_scaled * gradC * gradC.transpose();
    return hess;
}

} // namespace Energy