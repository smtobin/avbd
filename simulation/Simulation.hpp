#pragma once

#include "common/common.hpp"

#include "simulation/SimulationLogger.hpp"

#include "graphics/GraphicsScene.hpp"

#include "config/SimulationConfig.hpp"
#include "config/SimulationRenderConfig.hpp"

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

    void _timeStep();

    void _updateGraphics();

    protected:
    bool _setup;

    Real _time;
    Real _dt;
    Real _end_time;
    Real _g_accel;
    int _viewer_refresh_time_ms;

    Real _last_collision_check_time;


    int _solver_iters = 1;

    std::deque<std::function<void()>> _callback_queue;

    /** Responsible for logging various simulation quantities. */
    std::unique_ptr<SimulationLogger> _logger;

    // graphics
    Graphics::GraphicsScene _graphics_scene;

    Config::SimulationConfig _config;
};

} // namespace Simulation
