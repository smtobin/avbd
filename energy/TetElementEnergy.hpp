#pragma once

#include "energy/EnergyBase.hpp"
#include "common/Particle.hpp"
#include "common/mesh/TetMesh.hpp"

#include <array>

namespace Energy
{

class TetElementEnergy
{
public:
    TetElementEnergy(const TetMesh* mesh, int element_index, Real lambda, Real mu);
    
    /** Returns the current energy given the current state. */
    Real energy() const override;

    /** Computes the gradient of the energy with respect to a particular particle. */
    Vec3r gradient(int index) const override;

    /** Computes the Hessian of the energy function with respect to a particular particle. */
    Mat3r hessian(int index) const override; 

protected:
    /** Pointer to the mesh that owns this element */
    const TetMesh* _mesh;
    /** Element index in the tet mesh */
    int _element;

    /** Element material properties (Lame params) */
    Real _lambda;
    Real _mu;
};


} // namespace Energy