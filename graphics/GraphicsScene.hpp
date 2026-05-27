#pragma once

#include "common/mesh/ParticleMesh.hpp"


#include "config/SimulationRenderConfig.hpp"
#include "config/MeshRenderConfig.hpp"

#include "graphics/GraphicsObject.hpp"

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkOpenGLRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>

#include <vector>
#include <deque>
#include <atomic>
#include <memory>

namespace Sim
{
    class Simulation;
}

namespace Graphics
{

class GraphicsScene
{
    public:
    static void renderCallback(vtkObject* caller, long unsigned int event_id, void* client_data, void* call_data);
    void displayWindow() { _render_window->Render(); }
    void interactorStart() { _interactor->Start(); }

    explicit GraphicsScene();
    explicit GraphicsScene(const Config::SimulationRenderConfig& sim_render_config);

    void setup(Sim::Simulation* sim=nullptr);

    void update();

    void addObject(const SimObject::TetMeshObject* mesh_obj, const Config::MeshRenderConfig& render_config);
    void addObject(const SimObject::RigidSphere* /* sphere */, const Config::ObjectRenderConfig& /* render_config */) {}

    Vec3r cameraPosition() const;
    Vec3r cameraViewDirection() const;
    Vec3r cameraUpDirection() const;
    Vec2r worldCoordinatesToPixelCoordinates(const Vec3r& world) const;

private:
    vtkSmartPointer<vtkOpenGLRenderer> _renderer;
    vtkSmartPointer<vtkRenderWindow> _render_window;
    vtkSmartPointer<vtkRenderWindowInteractor> _interactor;

    std::vector<std::unique_ptr<GraphicsObject>> _graphics_objects;

    /** Stores meshes that are used as additional visualization representations of other objects.
     * I.e. objects may have multiple meshes associated with them for visualization purposes
     * 
     * Store in a deque so that pointers to meshes are valid upon resize
     */
    // std::deque<Mesh> _graphics_meshes;

    std::atomic<bool> _should_render;

    Config::SimulationRenderConfig _render_config;

};

} // namespace Graphics