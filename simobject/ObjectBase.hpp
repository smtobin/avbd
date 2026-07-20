#pragma once

#include "common/common.hpp"
#include "config/ObjectConfig.hpp"

#include <functional>
#include <atomic>

namespace SimObject
{

/** Base class for all objects in the sim */
class Object_Base
{
protected:
    /** Global object ID counter */
    inline static std::atomic<unsigned> _next_id{0};

    /** Pointer to the simulation context, which contains all simulation information */
    Sim::SimulationContext* _ctx;

    /** Object name */
    std::string _name;

    /** Object unique ID */
    unsigned _id;


public:
    Object_Base(Sim::SimulationContext* ctx, const Config::ObjectConfig& config)
        : _ctx(ctx)
        , _name(config.name())
        , _id(_next_id.fetch_add(1))
    {}

    Object_Base(const Object_Base&) = delete;
    Object_Base& operator=(const Object_Base&) = delete;

    Object_Base(Object_Base&&) = default;
    Object_Base& operator=(Object_Base&&) = default;

    virtual ~Object_Base() = default;

    /** Name of the object */
    const std::string& name() const { return _name; }

    /** Responsible for the initial setup of the object. */
    virtual void setup() = 0;

    /** Object ID */
    unsigned id() const { return _id; }

};

} // namespace SimObject