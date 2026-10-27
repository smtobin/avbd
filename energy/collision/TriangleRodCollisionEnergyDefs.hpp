#pragma once

#include "energy/collision/CollisionConstraintEnergyDefs.hpp"

namespace Energy
{

/** Required data for triangle-rod collision energy */
struct TriangleRodCollisionEnergyInfo : CollisionConstraintEnergyInfo
{
    Vec5u particle_indices;     // particle indices - first 3 indices are the triangle vertices, last 2 indices make up the rod segment
    Vec3r barys;    // barycentric coordinates of the contact point on the face
    Real s;  // location of contact point along rod segment (in the interval [0,1])
    Vec3r cp_rod_local; // local contact point in the rod interpolated frame
};

///////////////////////////////////////////////////////////////////////////

/** A non-owning view into the triangle-rod collision energy storage */
struct TriangleRodCollisionEnergyPoolView : CollisionConstraintEnergyPoolView<TriangleRodCollisionEnergyInfo>
{
};

} // namespace Energy