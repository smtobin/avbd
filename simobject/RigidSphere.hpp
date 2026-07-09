#include "simobject/ObjectBase.hpp"

#include "config/RigidSphereConfig.hpp"

namespace SimObject
{

class RigidSphere : public Object_Base
{
public: 
    RigidSphere(Sim::SimulationContext* ctx, const Config::RigidSphereConfig& config)
        : Object_Base(ctx, config)
    {}

    virtual void setup() override {}
    // virtual void for_each_particle(std::function<void(Particle*)>) override {}
    // virtual void for_each_particle(std::function<void(const Particle*)>) const override {}

    /** Provides a way to iterate through all energies owned by the object. */
    // virtual void for_each_energy(std::function<void(Energy::Energy_Base*)>) override {}
    // virtual void for_each_energy(std::function<void(const Energy::Energy_Base*)>) const override {}

    /** Provides a way to iterate through all QUADRATIC energies owned by the object. */
    // virtual void for_each_quadratic_energy(std::function<void(QuadraticEnergy*)>) override {}
    // virtual void for_each_quadratic_energy(std::function<void(QuadraticEnergy*)>) const override {}
};

} // namespace SimObject