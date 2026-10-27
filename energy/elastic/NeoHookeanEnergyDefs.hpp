#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"

namespace Energy
{

/** Required info for each NeoHookean energy. */
struct NeoHookeanEnergyInfo
{
    Vec4u particle_indices;     // indices of particles in the element
    Real lambda;                // first Lame parameter
    Real mu;                    // second Lame parameter
    Real kd;                    // strain-rate damping coefficient
    Mat3r Q;                    // rest state matrices (F = XQ^T)
    Real rest_volume;           // rest-state volume of the elements
    Mat3r E_prev;               // previous Green strain (used for strain-rate damping)
};

/////////////////////////////////////////////////////////////////////////////////////////////


/** A non-owning view into NeoHookean energy storage */
struct NeoHookeanEnergyPoolView : TombstonePoolView
{
    using PoolType = NeoHookeanEnergyPool;

    NeoHookeanEnergyInfo* data;
};


} // namespace Energy