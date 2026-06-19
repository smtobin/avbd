#include "simulation/Simulation.hpp"
#include "common/Math.hpp"

#include <chrono>
#include <thread>
#include <filesystem>

namespace Sim
{

Simulation::Simulation()
    : _setup(false),
      _time(0),
    //   , _dt(1e-3), _end_time(10),
    //   _g_accel(9.81), _viewer_refresh_time_ms(1000.0/30.0),
      _graphics_scene()
{
    _last_collision_check_time = std::numeric_limits<Real>::lowest();
}

Simulation::Simulation(const Config::SimulationConfig& sim_config)
    : _setup(false)
    , _time(0)
    // , _dt(sim_config.timeStep())
    // , _end_time(sim_config.endTime())
    // , _g_accel(sim_config.gAccel())
    // , _viewer_refresh_time_ms(1000.0/30.0)
    , _ctx(1000, 5000)
    , _solver(&_ctx, sim_config.solverIters(), sim_config.iterAcceleration())
    , _graphics_scene(sim_config.renderConfig())
    , _config(sim_config)
{
    _ctx.params.dt = sim_config.timeStep();
    _ctx.params.end_time = sim_config.endTime();
    _ctx.params.g_accel = sim_config.gAccel();
    _ctx.params.viewer_refresh_time_ms = 1000.0/30.0;   // 30 fps

    _last_collision_check_time = std::numeric_limits<Real>::lowest();
}

void Simulation::setup()
{
    _setup = true;

    // setup the graphics scene
    _graphics_scene.setup(this);

    if (_config.groundPlane())
    {
        _graphics_scene.addGroundPlane();
    }

    // set up the logger (when applicable)
    if (_config.logging())
    {
        // get datetime string
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d_%H:%M:%S") << ".txt";
        std::string filename = ss.str();

        std::filesystem::path output_dir(_config.loggingOutputDir());
        std::filesystem::path filepath = output_dir / filename;

        // create the logger
        _logger = std::make_unique<SimulationLogger>(filepath.string(), _config.loggingInterval());

        _logger->addOutput("time", &_time);
    }

    // create objects
    const ObjectConfigs_Container& obj_configs = _config.objectConfigs();
    obj_configs.for_each_element([&](const auto& obj_config){
        _addObjectFromConfig(obj_config); 
    });

    // after creating the objects, build the adjacency structure
    _ctx.adjacency.buildAdjacency(_ctx.particles, _ctx.energies);
    
    
}

void Simulation::update()
{
    // we assume that other derived Simulation classes have already added their logged quantities
    // so we can start logging now (which will print the header and prevent us from adding new logged quantities)
    if (_logger)
        _logger->startLogging();

    auto wall_time_start = std::chrono::steady_clock::now();
    auto last_redraw = std::chrono::steady_clock::now();

    while (_time < _ctx.params.end_time)
    {
        // run any callbacks that have been queued
        for (; !_callback_queue.empty(); _callback_queue.pop_front())
        {
            _callback_queue.front()();
        }
        
        // the elapsed seconds in wall time since the simulation has started
        Real wall_time_elapsed_s = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - wall_time_start).count() / 1000000000.0;
        
        // if the simulation is ahead of the current elapsed wall time, stall
        if (_config.simMode() == Config::SimulationMode::VISUALIZATION && _time > wall_time_elapsed_s)
        {
            continue;
        }

        // const XPBDCollisionConstraints_Container& new_collision_constraints = _collision_scene.detectCollisions();
        // new_collision_constraints.for_each_element([&](const auto& collision_constraint) {
        //     using ConstraintType = std::remove_cv_t<std::remove_reference_t<decltype(collision_constraint)>>;

        //     auto& constraint_vec = _constraints.template get<ConstraintType>();
        //     constraint_vec.emplace_back(collision_constraint);
        //     ConstVectorHandle<ConstraintType> constraint_ref(&constraint_vec, constraint_vec.size()-1);
        //     _solver.addConstraint(constraint_ref);
        // });

        _timeStep();

        // the time in ms since the viewer was last redrawn
        auto time_since_last_redraw_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - last_redraw).count();
        // we want ~30 fps, so update the viewer every 33 ms
        if (time_since_last_redraw_ms > _ctx.params.viewer_refresh_time_ms)
        {
            // std::cout << _haptic_device_manager->getPosition() << std::endl;
            _updateGraphics();

            last_redraw = std::chrono::steady_clock::now();
        }
    }

    if (_logger)
    {
        _logger->stopLogging();
    }

    auto wall_time_end = std::chrono::steady_clock::now();
    std::cout << "Simulation " << _ctx.params.end_time << " seconds took " << std::chrono::duration_cast<std::chrono::milliseconds>(wall_time_end - wall_time_start).count() << " ms" << std::endl;
}

int Simulation::run()
{
    // setup if we haven't already
    if (!_setup)
        setup();
    

    // start update thread
    _graphics_scene.displayWindow();

    std::thread update_thread;
    if (_config.simMode() != Config::SimulationMode::FRAME_BY_FRAME)
        update_thread = std::thread(&Simulation::update, this);
    
    _graphics_scene.interactorStart();
}

void Simulation::notifyKeyPressed(const std::string& /*key*/)
{
    if (_config.simMode() == Config::SimulationMode::FRAME_BY_FRAME)
    {
        _timeStep();
        _updateGraphics();
    }
}

void Simulation::notifyKeyReleased(const std::string& /*key*/)
{
    
}

void Simulation::notifyMouseMoved(double /*mx*/, double /*my*/)
{
    
}

void Simulation::notifyLeftMouseButtonPressed()
{

}

void Simulation::notifyLeftMouseButtonReleased()
{

}

void Simulation::_timeStep()
{
    std::cout << "t=" << _time << std::endl;

    // let the solver do the iterations
    _solver.solve(_ctx.params.dt);

    // log quantities
    if (_logger)
    {
        _logger->logToFile(_time);
    }

    _time += _ctx.params.dt;
}

void Simulation::_updateGraphics()
{
    _graphics_scene.update();
}

} // namespace Sim