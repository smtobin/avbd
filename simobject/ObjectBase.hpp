#pragma once

#include "common/common.hpp"
#include "common/Particle.hpp"
#include "config/ObjectConfig.hpp"

#include <functional>

namespace Sim
{
    class SimulationContext;
}

namespace SimObject
{

class Object_Base
{
public:
    Object_Base(Sim::SimulationContext* ctx, const Config::ObjectConfig& config)
        : _ctx(ctx), _name(config.name())
    {}

    virtual ~Object_Base() = default;

    const std::string& name() const { return _name; }

    /** Responsible for the initial setup of the object. */
    virtual void setup() = 0;

    /** Provides a way to iterate through the particles owned by the object. */
    // virtual void for_each_particle(std::function<void(Particle*)> func) = 0;
    // virtual void for_each_particle(std::function<void(const Particle*)> func) const = 0;

    /** Provides a way to iterate through all energies owned by the object. */
    // virtual void for_each_energy(std::function<void(Energy::Energy_Base*)> func) = 0;
    // virtual void for_each_energy(std::function<void(const Energy::Energy_Base*)> func) const = 0;

    /** Provides a way to iterate through all QUADRATIC energies owned by the object. */
    // virtual void for_each_quadratic_energy(std::function<void(QuadraticEnergy*)> func) = 0;
    // virtual void for_each_quadratic_energy(std::function<void(QuadraticEnergy*)> func) const = 0;

protected:
    Sim::SimulationContext* _ctx;

    std::string _name;
};

} // namespace SimObject