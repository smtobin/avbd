#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"

namespace Energy
{

/** Required data for Cosserat rod energy */
struct CosseratRodEnergyInfo
{
    Vec2u particle_indices;     // the indices of the oriented particles in the rod element
    Vec6r stiffness;            // the diagonal of the body-frame stiffness matrix (GA, GA, EA, EI, EI, GJ)     /** TODO: (08/25/26) Move this to rod collision geometry with radius? */
    Real rest_length;           // initial length of the element
    Vec3r precurvature;         // initial curvature of the element
};

/////////////////////////////////////////////////////////////////////////////////////////

/** Non-owning view into Cosserat energy pool */
struct CosseratRodEnergyPoolView
{
    CosseratRodEnergyInfo* data;
};

} // namespace Energy