#pragma once

#include "config/ObjectConfig.hpp"

namespace Config
{

class RodConfig : public ObjectConfig
{
public:
    using SimObjectType = SimObject::Rod;

    explicit RodConfig()
        : ObjectConfig()
    {}

    explicit RodConfig(const YAML::Node& node)
        : ObjectConfig(node)
    {
        _extractParameter("base-fixed", node, _base_fixed);
        _extractParameter("tip-fixed", node, _tip_fixed);

        _extractParameter("length", node, _length);
        _extractParameter("diameter", node, _diameter);
        _extractParameter("elements", node, _elements);

        _extractParameter("density", node, _density);
        _extractParameter("E", node, _E);
        _extractParameter("nu", node, _nu);

        _extractParameter("curvature", node, _curvature);        
    }

    bool baseFixed() const { return _base_fixed.value; }
    bool tipFixed() const { return _tip_fixed.value; }

    Real length() const { return _length.value; }
    Real diameter() const { return _diameter.value; }
    int elements() const { return _elements.value; }

    Real density() const { return _density.value; }
    Real E() const { return _E.value; }
    Real nu() const { return _nu.value; }

    Vec3r curvature() const { return _curvature.value; }

protected:
    ConfigParameter<bool> _base_fixed = ConfigParameter<bool>(true);
    ConfigParameter<bool> _tip_fixed = ConfigParameter<bool>(false);

    ConfigParameter<Real> _length = ConfigParameter<Real>(1.0);
    ConfigParameter<Real> _diameter = ConfigParameter<Real>(0.1);
    ConfigParameter<int> _elements = ConfigParameter<int>(20);
    
    ConfigParameter<Real> _density = ConfigParameter<Real>(1000);
    ConfigParameter<Real> _E = ConfigParameter<Real>(3e6);
    ConfigParameter<Real> _nu = ConfigParameter<Real>(0.45);

    ConfigParameter<Vec3r> _curvature = ConfigParameter<Vec3r>(Vec3r(0,0,0));
};
    
} // namespace Config