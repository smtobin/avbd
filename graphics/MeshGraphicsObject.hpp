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

    virtual void update() override;

private:
    /** The mesh geometry */
    const ParticleMesh* _mesh;

    /** Poly data for the mesh */
    vtkSmartPointer<vtkPolyData> _vtk_poly_data;

    /** Normals filter */
    vtkSmartPointer<vtkPolyDataNormals> _vtk_normals;
};

} // namespace Graphics