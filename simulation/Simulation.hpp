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
#include "simobject/Rod.hpp"
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

        unsigned prim_idx = _ctx.collision_pool.addObject(*new_obj_ptr);
        new_obj_ptr->setCollisionPrimitiveIndex(prim_idx);

        // add static ground-object collision constraints
        _addGroundCollisionConstraintsForObject(*new_obj_ptr);


        _graphics_scene.addObject(new_obj_ptr.get(), obj_config.renderConfig());

        _objects.push_back(std::move(new_obj_ptr));
    }
       
    void _addObjectFromConfig(const Config::TetMeshObjectConfig& obj_config)
    {
        // using ObjPtrType = std::unique_ptr<typename ConfigType::SimObjectType>;
        using ObjType = SimObject::TetMeshObject;
        using ObjPtrType = std::unique_ptr<ObjType>;
        ObjPtrType new_obj_ptr = std::make_unique<ObjType>(&_ctx, obj_config);
        new_obj_ptr->setup();

        _ctx.collision_pool.addObject(*new_obj_ptr);
        _addGroundCollisionConstraintsForObject(*new_obj_ptr);

        _graphics_scene.addObject(new_obj_ptr.get(), obj_config.meshRenderConfig());

        _objects.push_back(std::move(new_obj_ptr));
    }

    void _addObjectFromConfig(const Config::RodConfig& obj_config)
    {
        using ObjType = SimObject::Rod;
        using ObjPtrType = std::unique_ptr<ObjType>;
        ObjPtrType new_obj_ptr = std::make_unique<ObjType>(&_ctx, obj_config);
        new_obj_ptr->setup();

        unsigned sdf_slot = _ctx.collision_pool.addObject(*new_obj_ptr);
        new_obj_ptr->setSdfIndex(sdf_slot);

        _addGroundCollisionConstraintsForObject(*new_obj_ptr);

        _graphics_scene.addObject(new_obj_ptr.get(), obj_config.renderConfig());

        _objects.push_back(std::move(new_obj_ptr));
    }

    /** Helpers for adding ground collision constraints for each type of object */
    void _addGroundCollisionConstraintsForObject(const SimObject::RigidSphere& sphere)
    {
        /** TODO: (08/10/26) set coefficients of friction based on material properties */
        // initial energy stiffness should depend on the particle inertia
        Real k_start = _ctx.particles.masses[sphere.com()] / (_ctx.params.dt * _ctx.params.dt);
        unsigned sdf_index = _ctx.collision_pool.particle_indices[sphere.collisionPrimitiveIndex()][0];
        _ctx.energies.rigid_body_ground_collision.addEnergy(sphere.com(), sdf_index, Vec3r(0, -sphere.radius(), 0), k_start, 0.5, 0.2);
    }

    void _addGroundCollisionConstraintsForObject(const SimObject::TetMeshObject& tet_mesh_obj)
    {
        // add ground collision constraints for each particle in the mesh
        /** TODO: (08/05/26) set coefficients of friction based on material properties */
        for (auto& v_idx : tet_mesh_obj.mesh().vertices())
        {
            Real k_start = _ctx.particles.masses[v_idx] / (_ctx.params.dt * _ctx.params.dt);
            std::cout << "k start: " << k_start << std::endl;
            _ctx.energies.ground_collision.addEnergy(v_idx, k_start, 0.4, 0.2);
        }
    }

    void _addGroundCollisionConstraintsForObject(const SimObject::Rod& rod)
    {
        /** TODO: (08/24/26) Add ground collision constraints for rods. */
        for (unsigned n_idx : rod.nodes())
        {
            Real k_start = _ctx.particles.masses[n_idx] / (_ctx.params.dt * _ctx.params.dt);
            _ctx.energies.rigid_body_ground_collision.addEnergy(n_idx, rod.sdfIndex(), Vec3r(0, -rod.radius(), 0), k_start, 0.5, 0.2);
        }
        
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
