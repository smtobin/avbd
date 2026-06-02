#pragma once

#include "common/common.hpp"
#include "energy/EnergyBase.hpp"

namespace Energy
{

/** Base class for energy expressions of the form E = 1/2 * C(x)^T * K * C(x)  */
template<int ConstraintDim>
class QuadraticEnergy : public Energy_Base
{
public:
    using ConstraintVecType = Eigen::Vector<Real, ConstraintDim>;

    QuadraticEnergy(const ConstraintVecType& K_max, const ConstraintVecType& K_init = 1e2*ConstraintVecType::Ones())
        : _k_cur(K_init), _k_max(K_max), _k_start(K_init)
    {
        _last_C = ConstraintVecType::Zero();
    }

    virtual ~QuadraticEnergy() = default;

    virtual Real energy() const override
    {
        ConstraintVecType C = evaluateConstraint();
        return (0.5 * C.transpose() * _k_cur.asDiagonal() * C).value();
    }

    /** Resets the constraint stiffness to starting k
     * Should be called before the first iteration at each time step.
     */
    virtual void reset() override 
    { 
        // _k_cur = _k_start.array().max(STIFFNESS_GAMMA*_k_cur.array());
     }

    /** Update the constraint stiffness */
    virtual void updateAfterIteration() override
    {
        // ConstraintVecType C = evaluateConstraint();
        // _k_cur += STIFFNESS_BETA * C.cwiseAbs();
        // // _k_cur *= STIFFNESS_BETA;
        // _k_cur = _k_cur.array().min(_k_max.array());

        // _last_C = C;

        // std::cout << "New k: " << _k_cur.value() << std::endl;
    }

    /** Evaluate C(x) */
    virtual ConstraintVecType evaluateConstraint() const = 0;

protected:
    /** Current stiffness */
    ConstraintVecType _k_cur;

    /** Max stiffness */
    ConstraintVecType _k_max;

    /** Initial stiffness at each iteration */
    ConstraintVecType _k_start;

    /** Store the last constraint evaluation */
    ConstraintVecType _last_C;

};

} // namespace Energy