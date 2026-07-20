#pragma once

#include "simobject/rigid/RigidObjectBase.hpp"

#include "config/RigidSphereConfig.hpp"

namespace SimObject
{

class RigidSphere : public RigidObject_Base
{
private:
    Real _radius;

public: 
    RigidSphere(Sim::SimulationContext* ctx, const Config::RigidSphereConfig& config);

    virtual void setup() override {}

    Real radius() const { return _radius; }
};

} // namespace SimObject