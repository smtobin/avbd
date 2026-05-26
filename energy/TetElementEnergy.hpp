#pragma once

#include "energy/EnergyBase.hpp"

#include "common/Particle.hpp"

#include <array>

namespace Energy
{

class TetElementEnergy
{
public:
    TetElementEnergy(const std::array<const Particle*, 4>& element_particles);
    
    /** Returns the current energy given the current state. */
    Real energy() const override;

    /** Computes the gradient of the energy with respect to a particular particle. */
    Vec3r gradient(int index) const override;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    Mat3r hessian(int index) const override; 

protected:
    std::array<const Particle*, 4> _particles;
};


} // namespace Energy