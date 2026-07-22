#include "common/common.hpp"
#include "common/Algorithm.hpp"
#include "collision/LBVH.hpp"
#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"

#include "simobject/TetMeshObject.hpp"
#include "simobject/rigid/RigidSphere.hpp"


#include "simulation/SimulationContext.hpp"

#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkOpenGLRenderer.h>

#include <gmsh.h>

#include <bitset>

void testRadixTree()
{
    int capacity = 1000;
    Collision::CollisionPrimitivePool col_pool(capacity, 0);

    Collision::LBVH bvh;
    /** Test radix tree construction using the example from Karras 2012 */
    std::array<uint64_t, 8> morton_codes = {
        0b11001ULL,
        0b00100ULL,
        0b00101ULL,
        0b00010ULL,
        0b11000ULL,
        0b00001ULL,
        0b10011ULL,
        0b11110ULL
    };

    for (unsigned i = 0; i < morton_codes.size(); i++)
    {
        unsigned new_slot = col_pool.allocSlot();
        std::cout << "New slot: " << new_slot << std::endl;
        col_pool.morton_code[new_slot] = morton_codes[i];
    }

    Algorithm::radixSort(col_pool.morton_code, col_pool.sorted_order, col_pool.totalSize());
    std::cout << "Sorted order:" << std::endl;
    for (unsigned i = 0; i < col_pool.totalSize(); i++)
    {
        std::cout << " " << i << ": " << col_pool.sorted_order[i] << std::endl;
    }
    Collision::LBVHBuilder::constructTree(col_pool, bvh);
    std::cout << "Root: " << bvh.root << std::endl;
    for (unsigned i = 0; i < bvh.parent.size(); i++)
    {
        std::cout << "Node " << i << ": Left=" << bvh.left[i] << " Right=" << bvh.right[i] << " Parent=" << bvh.parent[i] << " Leaf count=" << bvh.leaf_count[i] << std::endl;
    }
    bvh.printTreeWithInfo(col_pool);
}

void visualizeBVH(const Collision::LBVH& lbvh)
{
    // find max depth
    std::vector<unsigned> depths = lbvh.nodeDepths();

    unsigned max_depth = *std::max_element(depths.begin(), depths.end());

    std::cout << "Max depth: " << max_depth << std::endl;

    // build geometry for each level
    std::vector<vtkSmartPointer<vtkAppendPolyData>> appenders(max_depth + 1);

    for (auto& a : appenders)
        a = vtkSmartPointer<vtkAppendPolyData>::New();

    for (unsigned i = 0; i < depths.size(); ++i)
    {
        auto cube = vtkSmartPointer<vtkCubeSource>::New();

        cube->SetBounds(
            lbvh.min_x[i], lbvh.max_x[i],
            lbvh.min_y[i], lbvh.max_y[i],
            lbvh.min_z[i], lbvh.max_z[i]);

        cube->Update();

        // Vec3r min_( lbvh.min_x[i],
        //     lbvh.min_y[i],
        //     lbvh.min_z[i]);
        // Vec3r max_(lbvh.max_x[i],
        //     lbvh.max_y[i],
        //     lbvh.max_z[i]);
        // std::cout << "Depth " << depths[i] << " AABB " << min_.transpose() << " to " << max_.transpose() << std::endl;

        appenders[depths[i]]->AddInputData(cube->GetOutput());
    }

    for (auto& a : appenders)
        a->Update();

    // create render window and renderer
    vtkNew<vtkOpenGLRenderer> renderer;
    renderer->SetBackground(0.1, 0.1, 0.15);

    vtkNew<vtkRenderWindow> render_window;
    render_window->AddRenderer(renderer);
    render_window->SetSize(1200, 900);

    vtkNew<vtkRenderWindowInteractor> interactor;
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    interactor->SetInteractorStyle(style);
    interactor->SetRenderWindow(render_window);

    // render each level separately
    vtkNew<vtkNamedColors> colors;

    std::vector<std::array<double,3>> palette =
    {
        {1.0,0.0,0.0},
        {1.0,0.5,0.0},
        {1.0,1.0,0.0},
        {0.0,1.0,0.0},
        {0.0,1.0,1.0},
        {0.0,0.4,1.0},
        {0.7,0.0,1.0}
    };

    for (unsigned d = 0; d <= max_depth; ++d)
    {
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(appenders[d]->GetOutputPort());

        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);

        auto c = palette[d % palette.size()];

        actor->GetProperty()->SetColor(c[0], c[1], c[2]);

        // actor->GetProperty()->SetOpacity(0.15);
        actor->GetProperty()->SetRepresentationToWireframe();
        actor->GetProperty()->SetLineWidth(2.0);
        actor->GetProperty()->SetOpacity(1.0);

        renderer->AddActor(actor);
    }

    renderer->ResetCamera();
    render_window->Render();
    interactor->Initialize();
    interactor->Start();
}

void testFewTrianglesBVH()
{
    Collision::CollisionPrimitivePool col_pool(1000, 0);
    std::vector<std::array<Vec3r, 3>> triangle_vertices = {
        {
            Vec3r(0.0, 0.0, 0.0),
            Vec3r(1.0, 0.0, 0.5),
            Vec3r(0.0, 1.0, 0.0)
        },
        {
            Vec3r(1.0, 1.0, 1.0),
            Vec3r(1.2, 1.0, 1.2),
            Vec3r(1.4, 1.2, 1.2)
        },
        {
            Vec3r(0.6, 0.8, 0.8),
            Vec3r(0.4, 0.6, 0.6),
            Vec3r(0.8, 0.4, 0.6)
        },
        {
            Vec3r(0.7, 1.1, 1.3),
            Vec3r(1.1, 1.2, 1.4),
            Vec3r(0.6, 0.7, 0.7)
        }
    };

    std::vector<std::array<unsigned, 3>> triangles;

    Sim::SimulationContext ctx;
    for (const auto& tri : triangle_vertices)
    {
        std::array<unsigned, 3> new_tri;
        for (unsigned v_idx = 0; v_idx < 3; v_idx++)
        {
            unsigned idx = ctx.particles.addParticle(tri[v_idx], 0);
            new_tri[v_idx] = idx;
        }
        triangles.push_back(new_tri);

        unsigned c_idx = col_pool.allocSlot();
        col_pool.type[c_idx] = Collision::PrimitiveType::Triangle;
        col_pool.particle_indices[c_idx][0] = new_tri[0];
        col_pool.particle_indices[c_idx][1] = new_tri[1];
        col_pool.particle_indices[c_idx][2] = new_tri[2];
        col_pool.num_particles[c_idx] = 3;
    }

    Collision::LBVH lbvh;
    Collision::LBVHBuilder::buildBVH(ctx.particles, col_pool, lbvh);

    std::vector<std::pair<unsigned, unsigned>> collision_pairs;
    Collision::LBVHTraversal::traverseSelfIterative(lbvh, lbvh.root, collision_pairs);
    for (const auto& collision_pair : collision_pairs)
    {
        std::cout << "Potential collision between nodes " << collision_pair.first << " and " << collision_pair.second << std::endl;
        std::cout << "BVH Node A - Prim type: " << static_cast<unsigned>(col_pool.type[collision_pair.first]) << "  Particles: ";
        for (unsigned k = 0; k < col_pool.num_particles[collision_pair.first]; k++)
        {
            std::cout << col_pool.particle_indices[collision_pair.first][k] << ", ";
        }
        std::cout << std::endl;
        std::cout << "BVH Node B - Prim type: " << static_cast<unsigned>(col_pool.type[collision_pair.second]) << "  Particles: ";
        for (unsigned k = 0; k < col_pool.num_particles[collision_pair.first]; k++)
        {
            std::cout << col_pool.particle_indices[collision_pair.second][k] << ", ";
        }
        std::cout << std::endl;
    }    

    visualizeBVH(lbvh);
    
}

void testSpheresAndMeshBVH()
{
    Sim::SimulationContext ctx;
    std::vector<Vec3r> sphere_locs = {
        Vec3r(1.0, 0.5, 0.2),
        Vec3r(1.6, 0.5, 0.4),
        Vec3r(1.2, 0.4, 0.3)
    };
    std::vector<Real> sphere_rads = {
        0.05, 0.02, 0.04
    };

    std::vector<SimObject::RigidSphere> spheres;
    for (unsigned i = 0; i < sphere_locs.size(); i++)
    {
        Config::RigidSphereConfig sphere_config(
            "sphere",
            sphere_locs[i],
            Vec3r::Zero(),
            Vec3r::Zero(),
            Vec3r::Zero(),
            true,
            0.5,
            0.2,
            1000,
            false,
            sphere_rads[i]
        );
        spheres.emplace_back(&ctx, sphere_config);
        ctx.collision_pool.addObject(spheres.back());
    }

    // Config::TetMeshObjectConfig mesh_config("../resource/cube2.msh");
    // SimObject::TetMeshObject mesh_obj(&ctx, mesh_config);
    // mesh_obj.setup();

    // ctx.collision_pool.addObject(mesh_obj);

    std::vector<std::array<Vec3r, 3>> triangle_vertices = {
        {
            Vec3r(0.0, 0.0, 0.0),
            Vec3r(1.0, 0.0, 0.5),
            Vec3r(0.0, 1.0, 0.0)
        }
        // ,
        // {
        //     Vec3r(1.0, 1.0, 1.0),
        //     Vec3r(1.2, 1.0, 1.2),
        //     Vec3r(1.4, 1.2, 1.2)
        // }
        // ,
        // {
        //     Vec3r(0.6, 0.8, 0.8),
        //     Vec3r(0.4, 0.6, 0.6),
        //     Vec3r(0.8, 0.4, 0.6)
        // },
        // {
        //     Vec3r(0.7, 1.1, 1.3),
        //     Vec3r(1.1, 1.2, 1.4),
        //     Vec3r(0.6, 0.7, 0.7)
        // }
    };

    std::vector<std::array<unsigned, 3>> triangles;

    for (const auto& tri : triangle_vertices)
    {
        std::array<unsigned, 3> new_tri;
        for (unsigned v_idx = 0; v_idx < 3; v_idx++)
        {
            unsigned idx = ctx.particles.addParticle(tri[v_idx], 0);
            new_tri[v_idx] = idx;
        }
        triangles.push_back(new_tri);

        unsigned c_idx = ctx.collision_pool.allocSlot();
        ctx.collision_pool.type[c_idx] = Collision::PrimitiveType::Triangle;
        ctx.collision_pool.particle_indices[c_idx][0] = new_tri[0];
        ctx.collision_pool.particle_indices[c_idx][1] = new_tri[1];
        ctx.collision_pool.particle_indices[c_idx][2] = new_tri[2];
        ctx.collision_pool.num_particles[c_idx] = 3;
    }

    Collision::LBVHBuilder::buildBVH(ctx.particles, ctx.collision_pool, ctx.lbvh);
    
    std::cout << "Morton codes: " << std::endl;
    for (unsigned i = 0; i < ctx.collision_pool.totalSize(); i++)
    {
        std::cout << " " << i << ": " << std::bitset<64>(ctx.collision_pool.morton_code[i]) << std::endl;
    }

    std::cout << "Sorted order:" << std::endl;
    for (unsigned i = 0; i < ctx.collision_pool.totalSize(); i++)
    {
        std::cout << " " << i << ": " << std::bitset<64>(ctx.collision_pool.morton_code[ctx.collision_pool.sorted_order[i]]) << std::endl;
    }

    ctx.lbvh.printTreeWithInfo(ctx.collision_pool);

    std::vector<std::pair<unsigned, unsigned>> collision_pairs;
    Collision::LBVHTraversal::traverseSelfIterative(ctx.lbvh, ctx.lbvh.root, collision_pairs);
    for (const auto& collision_pair : collision_pairs)
    {
        if (ctx.collision_pool.object_id[collision_pair.first] == ctx.collision_pool.object_id[collision_pair.second])
            continue;

        unsigned p1 = ctx.collision_pool.sorted_order[collision_pair.first];
        unsigned p2 = ctx.collision_pool.sorted_order[collision_pair.second];
        std::cout << "Potential collision between nodes " << collision_pair.first << " and " << collision_pair.second << std::endl;
        std::cout << "BVH Node A - Prim type: " << static_cast<unsigned>(ctx.collision_pool.type[p1]) << "  Particles: ";
        for (unsigned k = 0; k < ctx.collision_pool.num_particles[p1]; k++)
        {
            std::cout << ctx.collision_pool.particle_indices[p1][k] << ", ";
        }
        std::cout << std::endl;
        std::cout << "BVH Node B - Prim type: " << static_cast<unsigned>(ctx.collision_pool.type[p2]) << "  Particles: ";
        for (unsigned k = 0; k < ctx.collision_pool.num_particles[p2]; k++)
        {
            std::cout << ctx.collision_pool.particle_indices[p2][k] << ", ";
        }
        std::cout << std::endl;
    }  

    visualizeBVH(ctx.lbvh);
}

void testTetMeshBVH()
{
    Sim::SimulationContext ctx;

    Config::TetMeshObjectConfig mesh_config("../resource/cube2.msh");
    SimObject::TetMeshObject mesh_obj(&ctx, mesh_config);
    mesh_obj.setup();

    Collision::CollisionPrimitivePool col_pool(1000, 0);
    col_pool.addObject(mesh_obj);

    Collision::LBVH lbvh;
    Collision::LBVHBuilder::buildBVH(ctx.particles, col_pool, lbvh);

    visualizeBVH(lbvh);
}

int main()
{
    gmsh::initialize();

    // testRadixTree();
    testTetMeshBVH();
    testFewTrianglesBVH();
    testSpheresAndMeshBVH();
}