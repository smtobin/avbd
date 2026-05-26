#pragma once

#include "common/common.hpp"

namespace Energy
{

/** Base class for potential energy functions that affect particles in the system.
 * Defines a common interface that all energy functions abide by.
 */
class Energy_Base
{
public:
    Energy_Base() = default;
    virtual ~Energy_Base() = default;

    /** Returns the current energy given the current state. */
    Real energy() const = 0;

    /** Computes the gradient of the energy with respect to a particular particle. */
    Vec3r gradient(int index) const = 0;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    Mat3r hessian(int index) const = 0;

protected:
};

} // namespace Energy