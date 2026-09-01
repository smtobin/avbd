#include "graphics/GraphicsScene.hpp"
#include "graphics/CustomVTKInteractorStyle.hpp"

#include "simulation/Simulation.hpp"
#include "graphics/MeshGraphicsObject.hpp"
#include "graphics/RodGraphicsObject.hpp"
#include "graphics/SphereGraphicsObject.hpp"
#include "graphics/PlaneGraphicsObject.hpp"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <vtkHDRReader.h>
#include <vtkImageFlip.h>
#include <vtkImageReader2.h>
#include <vtkImageReader2Factory.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLight.h>
#include <vtkNamedColors.h>
#include <vtkOpenGLRenderer.h>
#include <vtkOpenGLTexture.h>
#include <vtkPBRIrradianceTexture.h>
#include <vtkPNGReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataTangents.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSkybox.h>
#include <vtkTexture.h>
#include <vtkTriangleFilter.h>

#include <vtkSequencePass.h>
#include <vtkShadowMapBakerPass.h>
#include <vtkShadowMapPass.h>
#include <vtkCameraPass.h>
#include <vtkRenderPassCollection.h>
#include <vtkRenderStepsPass.h>
#include <vtkToneMappingPass.h>
#include <vtkLightsPass.h>
#include <vtkOpaquePass.h>
#include <vtkTranslucentPass.h>
#include <vtkImageShiftScale.h>
#include <vtkDepthPeelingPass.h>

#include <vtkCoordinate.h>

namespace Graphics
{

GraphicsScene::GraphicsScene(const Config::SimulationRenderConfig& render_config)
    : _should_render(false), _render_config(render_config)
{
}

GraphicsScene::GraphicsScene()
    : _should_render(false), _render_config()
{
}

void GraphicsScene::renderCallback(vtkObject* /*caller*/, long unsigned int /*event_id*/, void* client_data, void* /*call_data*/)
{
    
    GraphicsScene* scene = static_cast<GraphicsScene*>(client_data);
    // std::cout << "Callback! t=" << simulation->_time << std::endl;
    if (scene->_should_render.exchange(false))
    {
        // std::cout << "Rendering... t=" << simulation->_time << std::endl;
        // set the clipping range every time
        vtkCamera* camera = scene->_renderer->GetActiveCamera();
        camera->SetClippingRange(0.01, 10000.0);

        scene->_render_window->Render();
    }
}


void GraphicsScene::setup(Sim::Simulation* sim)
{
    // create renderer for actors in the scene
    _renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
    
    // add all the rods that may have been added before we set up
    for (const auto& graphics_obj : _graphics_objects)
    {
        _renderer->AddActor(graphics_obj->actor());
    }
    

    _renderer->SetBackground(0.8, 0.8, 0.8);
    _renderer->SetAutomaticLightCreation(false);

    // set camera to look at (0,0,0) from its current position
    _renderer->ResetCamera();
    vtkCamera* camera = _renderer->GetActiveCamera();
    double focal_point[3];
    camera->GetFocalPoint(focal_point);

    // set focal point
    camera->SetFocalPoint(
        _render_config.cameraFocalPoint()[0],
        _render_config.cameraFocalPoint()[1],
        _render_config.cameraFocalPoint()[2]
    );

    // position camera
    camera->SetPosition(
        _render_config.cameraPosition()[0],
        _render_config.cameraPosition()[1],
        _render_config.cameraPosition()[2]
    );

    camera->SetViewUp(0, 1, 0);
    camera->SetClippingRange(0.01, 1000.0);

    camera->SetParallelProjection(_render_config.cameraOrthographic());

    //////////////////////////////////////////////////////////
    // Create HDR lighting (if specified in the config)
    /////////////////////////////////////////////////////////

    std::optional<std::string> hdr_filename = _render_config.hdrImageFilename();
    if (hdr_filename.has_value())
    {
        vtkNew<vtkTexture> hdr_texture;
        vtkNew<vtkHDRReader> reader;
        reader->SetFileName(hdr_filename.value().c_str());
        vtkNew<vtkImageShiftScale> scale;
        scale->SetInputConnection(reader->GetOutputPort());
        scale->SetScale(_render_config.hdrScaling()); // darker HDRI

        hdr_texture->SetInputConnection(scale->GetOutputPort());
        hdr_texture->SetColorModeToDirectScalars();
        hdr_texture->MipmapOn();
        hdr_texture->InterpolateOn();
        hdr_texture->SetInputConnection(scale->GetOutputPort());

        if (_render_config.createSkybox())
        {
            vtkNew<vtkSkybox> skybox;
            skybox->SetTexture(hdr_texture);
            skybox->SetFloorRight(0,0,1);
            skybox->SetProjection(vtkSkybox::Sphere);
            _renderer->AddActor(skybox);
        }

        _renderer->UseImageBasedLightingOn();
        _renderer->UseSphericalHarmonicsOff();
        _renderer->SetEnvironmentTexture(hdr_texture, false);
    }
    

    ////////////////////////////////////////////////////////
    // Add lights
    ////////////////////////////////////////////////////////

    vtkNew<vtkLight> light;
    light->SetLightTypeToSceneLight();
    light->SetPositional(true);
    light->SetPosition(0.0, 10.0, 0.0);
    light->SetFocalPoint(0, 0, 0);
    light->SetConeAngle(45);             
    light->SetColor(1.0, 1.0, 1.0);
    light->SetIntensity(1.0);
    light->SetShadowAttenuation(0.5);
    _renderer->AddLight(light);


    //////////////////////////////////////////////////////
    // Create the render window and interactor
    //////////////////////////////////////////////////////
    _render_window = vtkSmartPointer<vtkRenderWindow>::New();
    _render_window->AddRenderer(_renderer);
    _render_window->SetSize(600, 600);
    _render_window->SetWindowName("Rod Test");

    _interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    // vtkNew<vtkInteractorStyleTrackballCamera> style;
    vtkNew<CustomVTKInteractorStyle> style;
    style->AutoAdjustCameraClippingRangeOff();
    style->registerSimulation(sim);
    _interactor->SetInteractorStyle(style);
    _interactor->SetRenderWindow(_render_window);
    
    
    
    /////////////////////////////////////////////////////////
    // Create the rendering passes and settings
    ////////////////////////////////////////////////////////
    if (_render_config.fancyGraphics())
    {
        _render_window->SetAlphaBitPlanes(1);
        // _render_window->SetMultiSamples(0);

        _renderer->SetUseDepthPeeling(true);
        _renderer->SetMaximumNumberOfPeels(50);
        _renderer->SetOcclusionRatio(0.1);

        // Core passes
        vtkNew<vtkRenderPassCollection> passes;

        vtkNew<vtkLightsPass> lightsPass;
        vtkNew<vtkOpaquePass> opaquePass;
        vtkNew<vtkTranslucentPass> translucentPass;

        // Shadow map
        vtkNew<vtkShadowMapPass> shadowPass;
        shadowPass->GetShadowMapBakerPass()->SetResolution(4096);

        // IMPORTANT: depth peeling must wrap translucent pass
        vtkNew<vtkDepthPeelingPass> peelingPass;
        peelingPass->SetTranslucentPass(translucentPass);
        peelingPass->SetMaximumNumberOfPeels(50);
        peelingPass->SetOcclusionRatio(0.1);

        // Order matters
        passes->AddItem(shadowPass->GetShadowMapBakerPass());
        passes->AddItem(lightsPass);
        passes->AddItem(opaquePass);
        passes->AddItem(shadowPass);
        passes->AddItem(peelingPass);

        vtkNew<vtkSequencePass> seq;
        seq->SetPasses(passes);

        vtkNew<vtkCameraPass> cameraP;
        cameraP->SetDelegatePass(seq);

        _renderer->SetPass(cameraP);
    }

    vtkNew<vtkCallbackCommand> render_callback;
    render_callback->SetCallback(GraphicsScene::renderCallback);
    render_callback->SetClientData(this);
    _interactor->Initialize();
    _interactor->AddObserver(vtkCommand::TimerEvent, render_callback);
    _interactor->CreateRepeatingTimer(5);
}

void GraphicsScene::update()
{
    for (auto& obj : _graphics_objects)
    {
        obj->update();
    }

    _should_render.store(true);
}

void GraphicsScene::addObject(const SimObject::TetMeshObject* mesh_obj, const Config::MeshRenderConfig& render_config)
{
    std::unique_ptr<MeshGraphicsObject> mesh_go = std::make_unique<MeshGraphicsObject>(&mesh_obj->mesh(), render_config);
    _renderer->AddActor(mesh_go->actor());

    _graphics_objects.push_back(std::move(mesh_go));
}

void GraphicsScene::addObject(const SimObject::Rod* rod, const Config::ObjectRenderConfig& render_config)
{
    std::unique_ptr<RodGraphicsObject> rod_go = std::make_unique<RodGraphicsObject>(rod, render_config);
    _renderer->AddActor(rod_go->actor());

    _graphics_objects.push_back(std::move(rod_go));
}

void GraphicsScene::addObject(const SimObject::RigidSphere* sphere, const Config::ObjectRenderConfig& render_config)
{
    std::unique_ptr<SphereGraphicsObject> sphere_go = std::make_unique<SphereGraphicsObject>(sphere, render_config);
    _renderer->AddActor(sphere_go->actor());

    _graphics_objects.push_back(std::move(sphere_go));
}

void GraphicsScene::addGroundPlane()
{
    std::unique_ptr<PlaneGraphicsObject> plane_go = std::make_unique<PlaneGraphicsObject>(Vec3r::Zero(), Vec3r(0,1,0), 1000, 1000, Config::ObjectRenderConfig{});
    _renderer->AddActor(plane_go->actor());

    _graphics_objects.push_back(std::move(plane_go));
}

Vec3r GraphicsScene::cameraPosition() const
{
    double px, py, pz;
    _renderer->GetActiveCamera()->GetPosition(px, py, pz);
    return Vec3r(px, py, pz);
}

Vec3r GraphicsScene::cameraViewDirection() const
{
    double px, py, pz;
    _renderer->GetActiveCamera()->GetDirectionOfProjection(px, py, pz);
    return Vec3r(px, py, pz);

}

Vec3r GraphicsScene::cameraUpDirection() const
{
    double px, py, pz;
    _renderer->GetActiveCamera()->GetViewUp(px, py, pz);
    return Vec3r(px, py, pz);
}

Vec2r GraphicsScene::worldCoordinatesToPixelCoordinates(const Vec3r& world) const
{
    
    vtkNew<vtkCoordinate> coordinate;
    coordinate->SetCoordinateSystemToWorld();
    coordinate->SetValue(world[0], world[1], world[2]);

    int* display_coords = coordinate->GetComputedDisplayValue(_renderer);

    return Vec2r(display_coords[0], display_coords[1]);
}

} // namespace Graphics