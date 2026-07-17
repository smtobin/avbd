#pragma once

#include "config/RigidObjectConfig.hpp"

#include "simobject/ObjectBase.hpp"

namespace SimObject
{

class RigidObject_Base : public Object_Base
{
protected:
    /** Index in the oriented particle pool corresponding to the center of mass of the rigid body */
    unsigned _com;

public:
    RigidObject_Base(Sim::SimulationContext* ctx, const Config::RigidObjectConfig& config);

    unsigned com() const { return _com; }
    const Vec3r& position() const;
    const Quaternion& rotation() const;
    bool fixed() const;
};

} // namespace SimObject