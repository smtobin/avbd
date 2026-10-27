#pragma once

#include "common/common.hpp"
#include "common/TombstonePoolView.hpp"

namespace Energy
{

struct CollisionConstraintEnergyInfo
{
    Real k_start;                        // initial finite stiffness - this should scale depending on the inertia of the particle(s) affected
    Real k_n; Real k_t; Real k_b;        // the finite stiffnesses of the quadratic energy (one for each direction of the collision basis)
    Real lambda_n; Real lambda_t; Real lambda_b;   // the Lagrange multipliers enforcing the constraints (normal, tangent, binormal)
    Real C_n_prev;   // the constraint violation at the end of the previous time step (just for the normal direction)
    // (we do not use the constraint error smoothing for the tangent or binormal directions because the forces are naturally capped by the friction cone)

    /** The collision orthonormal frame */
    Vec3r normal;   // collision normal vector
    Vec3r tangent;  // collision tangent vector
    Vec3r binormal; // collision binormal vector

    /** Friction properties */
    bool use_static; // whether or not static friction is active
    // tracks whether the friction Lagrange multiplier was unclamped during the last iteration 
    // if, at the end of the time step, this is true and we were using kinetic friction, we switch to static
    bool unclamped_lambda_tb;   
    Real mu_s;  // static friction coefficient
    Real mu_k;  // kinetic friction coefficient
    
};

//////////////////////////////////////////////////////////////////////////////

/** A non-owning view into collision constraint energy pool. */
template <typename EnergyInfo>
struct CollisionConstraintEnergyPoolView : public TombstonePoolView
{
    EnergyInfo* data;
};


} // namespace Energy