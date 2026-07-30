#include "simobject/rigid/RigidSphere.hpp"

#include "simulation/SimulationContext.hpp"

namespace SimObject
{

RigidSphere::RigidSphere(Sim::SimulationContext* ctx, const Config::RigidSphereConfig& config)
    : RigidObject_Base(ctx, config)
    , _radius(config.radius())
{   
    // mass = 4/3 * pi * r^3 * density
    Real mass = 4 * M_PI / 3 * _radius * _radius * _radius * config.density();
    // moment of inertia = 2/5 * m * r^2
    Vec3r rot_inertia = 0.4 * mass * _radius * _radius * Vec3r::Ones();
    // rigid object base constructor creates particle, we need to set the inertia
    _ctx->particles.masses[_com] = mass;
    _ctx->particles.rotation_pool.rotational_inertias[_ctx->particles.rotation_idx[_com]] = rot_inertia;

}


} // namespace SimObject