#include "simobject/rigid/RigidObjectBase.hpp"

#include "simulation/SimulationContext.hpp"

#include "common/Math.hpp"

namespace SimObject
{

RigidObject_Base::RigidObject_Base(Sim::SimulationContext* ctx, const Config::RigidObjectConfig& config)
    : Object_Base(ctx, config)
{
    // create the COM particle
    // set the mass to be 0 initially - this is the responsibility of the derived rigid body classes to set
    Quaternion initial_rot = Math::QuaternionFromXYZEulerAngles(config.initialRotation());
    _com = _ctx->oriented_particles.addParticle(config.initialPosition(), initial_rot, 0, Vec3r::Zero());

    // set COM properties
    _ctx->oriented_particles.fixed[_com] = config.fixed();
}

const Vec3r& RigidObject_Base::position() const { return _ctx->oriented_particles.positions[_com]; }
const Quaternion& RigidObject_Base::rotation() const { return _ctx->oriented_particles.rotations[_com]; }
bool RigidObject_Base::fixed() const { return _ctx->oriented_particles.fixed[_com]; }
} // namespace SimObject