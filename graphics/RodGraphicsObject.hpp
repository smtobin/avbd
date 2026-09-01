#pragma once

#include "simobject/Rod.hpp"
#include "graphics/GraphicsObject.hpp"
#include "config/ObjectRenderConfig.hpp"

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>

#include <vtkActor.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>

namespace Graphics
{

class RodGraphicsObject : public GraphicsObject
{
public:
    explicit RodGraphicsObject(const SimObject::Rod* rod, const Config::ObjectRenderConfig& render_config);

    virtual void update() override;

private:
    void _generateInitialPolyData();
    void _updatePolyData();

    private:
    const SimObject::Rod* _rod;

    vtkSmartPointer<vtkPolyData> _vtk_poly_data;
    vtkSmartPointer<vtkPolyDataNormals> _vtk_poly_data_normals;
    vtkSmartPointer<vtkPolyDataMapper> _vtk_poly_data_mapper;

    Config::ObjectRenderConfig _render_config;

    /** Points (in the XY plane) for the cross-section of the rod */
    std::vector<Vec3r> _cross_section_points;

    /** Whether or not to color each rod element individually.
     * When enabled, a fixed number of samples per element will be used.
     */
    bool _color_elements;

    /** Number of points to sample for each element */
    int _sample_points_per_element = 3;

    /** Whether or not to visualize the centerline of the rod. */
    bool _draw_centerline;

    /** Whether or not to draw end caps on the rod. */
    bool _draw_end_caps;

    /** Number of cross sections in the end caps (when applicable) */
    int _cap_resolution = 6;

    /** Number of points in each circular cross section */
    int _tubular_resolution = 20;



};

} // namespace Graphics