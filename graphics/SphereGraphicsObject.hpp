#pragma once

#include "simobject/rigid/RigidSphere.hpp"
#include "config/ObjectRenderConfig.hpp"
#include "graphics/GraphicsObject.hpp"

#include <vtkSphereSource.h>
#include <vtkActor.h>
#include <vtkTransform.h>

namespace Graphics
{

class SphereGraphicsObject : public GraphicsObject
{
    public:
    explicit SphereGraphicsObject(const SimObject::RigidSphere* sphere, const Config::ObjectRenderConfig& render_config);

    virtual void update() override;

    const SimObject::RigidSphere* sphere() const { return _sphere; }

    private:
    const SimObject::RigidSphere* _sphere;

    vtkSmartPointer<vtkSphereSource> _sphere_source;    // sphere source algorithm that generates box mesh (keep a reference to it so we can potentially change it's size)
    vtkSmartPointer<vtkTransform> _vtk_transform;   // transform for the actor (updated according to position and orientation of sphere)
};

} // namespace Graphics