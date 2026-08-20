#pragma once

#include "simobject/ObjectBase.hpp"
#include "config/RodConfig.hpp"

namespace SimObject
{

class Rod : public Object_Base
{
public:
    Rod(Sim::SimulationContext* ctx, const Config::RodConfig& config);

    virtual void setup() override;

protected:
    std::vector<unsigned> _nodes;

    Real _radius;
    Real _length;

    /** Material properties */
    Real _E;
    Real _nu;
    Real _density;
};

} // namespace SimObject