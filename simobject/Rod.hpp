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
    unsigned _num_elements;

    Real _length;
    Vec3r _curvature;

    /** Material properties */
    Real _E;
    Real _nu;
    Real _density;
    Real _G;

    /** Cross-section properties */
    Real _radius;
    Real _area;
    Real _Ix;
    Real _Iy;
    Real _Iz;
};

} // namespace SimObject