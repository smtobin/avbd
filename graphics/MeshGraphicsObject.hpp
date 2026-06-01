#include "graphics/GraphicsObject.hpp"
#include "config/ObjectRenderConfig.hpp"
#include "config/MeshRenderConfig.hpp"

#include "common/mesh/ParticleMesh.hpp"

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkExtractEdges.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkTransform.h>

namespace Graphics
{

class MeshGraphicsObject : public GraphicsObject
{
public:
    explicit MeshGraphicsObject(const ParticleMesh* mesh, const Config::MeshRenderConfig& render_config);

    vtkSmartPointer<vtkActor> edgesActor() { return _edges_vtk_actor; }

    virtual void update() override;

private:
    /** The mesh geometry */
    const ParticleMesh* _mesh;

    /** Actor  for drawing edges of the mesh */
    vtkSmartPointer<vtkActor> _edges_vtk_actor;

    /** Poly data for the mesh */
    vtkSmartPointer<vtkPolyData> _vtk_poly_data;

    /** Transform for the mesh COM */
    vtkSmartPointer<vtkTransform> _vtk_transform;
};

} // namespace Graphics