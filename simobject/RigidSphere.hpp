#include "simobject/ObjectBase.hpp"

#include "config/RigidSphereConfig.hpp"

namespace SimObject
{

class RigidSphere : public Object_Base
{
public: 
    RigidSphere(const Config::RigidSphereConfig& config)
        : Object_Base(config)
    {}

    virtual void setup() override {}
    virtual void for_each_particle(std::function<void(Particle*)> func) override {}
    virtual void for_each_particle(std::function<void(const Particle*)> func) const override {}
};

} // namespace SimObject