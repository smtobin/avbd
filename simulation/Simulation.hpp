#pragma once

#include "common/common.hpp"

#include "simulation/SimulationContext.hpp"
// #include "simulation/SimulationLogger.hpp"
#include "simulation/SimulationExecutor.hpp"

#include "graphics/GraphicsScene.hpp"

#include "config/SimulationConfig.hpp"
#include "config/SimulationRenderConfig.hpp"

#include "collision/CollisionDetector.hpp"

#include "simobject/TetMeshObject.hpp"
#include "simobject/rigid/RigidSphere.hpp"

#include <vector>
#include <deque>
#include <unordered_map>


namespace Sim
{

class Simulation
{
    public:
    explicit Simulation();

    explicit Simulation(const Config::SimulationConfig& sim_config);

    virtual ~Simulation() = default;

    virtual void setup();

    virtual int run();

    virtual void update();


    // event handling
    virtual void notifyKeyPressed(const std::string& key);
    virtual void notifyKeyReleased(const std::string& key);
    virtual void notifyMouseMoved(double mx, double my);
    virtual void notifyLeftMouseButtonPressed();
    virtual void notifyLeftMouseButtonReleased();

    protected:

    template<typename ConfigType>        
    void _addObjectFromConfig(const ConfigType& obj_config)
    {
        // using ObjPtrType = std::unique_ptr<typename ConfigType::SimObjectType>;
        using ObjType = typename ConfigType::SimObjectType;
        using ObjPtrType = std::unique_ptr<ObjType>;

        ObjPtrType new_obj_ptr = std::make_unique<ObjType>(&_ctx, obj_config);
        new_obj_ptr->setup();

        _ctx.collision_pool.addObject(*new_obj_ptr);
        _graphics_scene.addObject(new_obj_ptr.get(), obj_config.renderConfig());

        _objects.push_back(std::move(new_obj_ptr));
    }
       
    void _addObjectFromConfig(const Config::TetMeshObjectConfig& obj_config)
    {
        // using ObjPtrType = std::unique_ptr<typename ConfigType::SimObjectType>;
        using ObjType = SimObject::TetMeshObject;
        using ObjPtrType = std::unique_ptr<ObjType>;
        // ObjType* new_obj_ptr = nullptr;
        // if constexpr (std::is_base_of_v<SimObject::ObjectGroup_Base, ObjType>)
        // {
        //     _object_groups.template push_back<ObjPtrType>(std::make_unique<ObjType>(obj_config));
        //     new_obj_ptr = _object_groups.template get<ObjPtrType>().back().get();
        //     new_obj_ptr->setup();

            

        //     // add the ObjectGroup's constraints to the solver
        //     _addConstraintsFromObject(new_obj_ptr, obj_config.projectorType());
        // }
        // else
        // {
            ObjPtrType new_obj_ptr = std::make_unique<ObjType>(&_ctx, obj_config);
            new_obj_ptr->setup();
            // _objects.template emplace_back<ObjPtrType>(std::make_unique<ObjType>(obj_config));
            // new_obj_ptr = _objects.template get<ObjPtrType>().back().get();
            // new_obj_ptr->setup();
        // }

        _ctx.collision_pool.addObject(*new_obj_ptr);
        _graphics_scene.addObject(new_obj_ptr.get(), obj_config.meshRenderConfig());

        // if the particles of this object should be logged, create logging outputs for them
        // if (obj_config.logParticles() && _logger)
        // {
        //     int particle_index = 0;
        //     new_obj_ptr->for_each_particle([&] (const Particle* particle) {
        //         const std::string var_name = new_obj_ptr->name() + "_particle" + std::to_string(particle_index++);
        //         _logger->addOutput(var_name, particle);
        //     });
        // }

        _objects.push_back(std::move(new_obj_ptr));
    }

    void _timeStep();

    void _updateGraphics();

    protected:
    bool _setup;

    Real _time;
    // Real _dt;
    // Real _end_time;
    // Real _g_accel;
    // int _viewer_refresh_time_ms;

    Real _last_collision_check_time;


    /** Number of solver iterations */
    // int _solver_iters = 1;

    /** Acceleration parameter "rho" for Chebyshev acceleration. VBD eqn (18) */
    // Real _iter_acceleration;

    std::deque<std::function<void()>> _callback_queue;

    /** The simulation context
     * Stores all particles, energies in the sim.
     */
    SimulationContext _ctx;

    /** Responsible for executing the simulation */
    SimulationExecutor _executor;

    /** Simulation objects
     * 
     * Mostly used for visualization and grouping the particles in the sim.
     */
    std::vector<std::unique_ptr<SimObject::Object_Base>> _objects;

    /** Responsible for logging various simulation quantities. */
    // std::unique_ptr<SimulationLogger> _logger;

    /** Responsible for visualization in the sim */
    Graphics::GraphicsScene _graphics_scene;

    /** Simulation params */
    Config::SimulationConfig _config;
};

} // namespace Sim
