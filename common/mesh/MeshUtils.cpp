#include "common/mesh/MeshUtils.hpp"

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags

#include <gmsh.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

ParticleMesh MeshUtils::loadSurfaceMeshFromFile(const std::string& filename)
{
    Assimp::Importer importer;
    importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS, 
        aiComponent_NORMALS | 
        aiComponent_TANGENTS_AND_BITANGENTS |
        aiComponent_COLORS |
        aiComponent_TEXCOORDS
    );

    // And have it read the given file with some example postprocessing
    // Usually - if speed is not the most important aspect for you - you'll
    // probably to request more postprocessing than we do in this example.
    const aiScene* scene = importer.ReadFile( filename,
        aiProcess_Triangulate            |
        aiProcess_JoinIdenticalVertices  |
        aiProcess_RemoveComponent        |
        aiProcess_SortByPType);

    // If the import failed, report it
    if (scene == nullptr)
    {
        std::cerr << "\tAssimp::Importer could not open " << filename << std::endl;
        std::cerr << "\tEnsure that the file is in a format that assimp can handle." << std::endl;
        assert(0);
    }

    const aiMesh* ai_mesh = scene->mMeshes[0];

    // Extract vertices
    std::vector<Vec3r> verts(ai_mesh->mNumVertices);
    for (unsigned i = 0; i < ai_mesh->mNumVertices; i++)
    {
        verts[i][0] = ai_mesh->mVertices[i].x;
        verts[i][1] = ai_mesh->mVertices[i].y;
        verts[i][2] = ai_mesh->mVertices[i].z;
    }

    // Extract faces
    std::vector<Vec3i> faces(ai_mesh->mNumFaces);
    for (unsigned i = 0; i < ai_mesh->mNumFaces; i++)
    {
        faces[i][0] = ai_mesh->mFaces[i].mIndices[0];
        faces[i][1] = ai_mesh->mFaces[i].mIndices[1];
        faces[i][2] = ai_mesh->mFaces[i].mIndices[2];
    }

    ParticleMesh mesh(verts, faces);
    return mesh;
}

ParticleTetMesh MeshUtils::loadTetMeshFromGmshFile(const std::string& filename)
{
    std::cout << "MeshUtils::loadTetMeshFromGmshFile - loading mesh data from " << filename << " as a ParticleMesh..." << std::endl;

    // ensure the file exists
    if (!std::filesystem::exists(filename))
    {
        std::cerr << "\t" << filename << " does not exist!" << std::endl;
        assert(0);
    }


    std::filesystem::path file_path(filename);

    if (file_path.extension() == ".obj" || file_path.extension() == ".stl")
    {
        convertToSTL(filename);
        convertSTLtoMSH(file_path.replace_extension(".stl"));
    }

    file_path = file_path.replace_extension(".msh");

    // ensure the file is a .msh file
    if (file_path.extension() != ".msh" && file_path.extension() != ".MSH")
    {
        std::cerr << "\t" << filename << " is not a .msh file!" << std::endl;
        assert(0);
    }

    gmsh::open(file_path.string());

    // Get all the elementary entities in the model, as a vector of (dimension,
    // tag) pairs:
    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities);

    std::vector<Vec3r> vertices;
    std::vector<Vec3i> faces;
    std::vector<Vec4i> elements;

    for(auto e : entities) {
        // Dimension and tag of the entity:
        int dim = e.first, tag = e.second;

        // Mesh data is made of `elements' (points, lines, triangles, ...), defined
        // by an ordered list of their `nodes'. Elements and nodes are identified by
        // `tags' as well (strictly positive identification numbers), and are stored
        // ("classified") in the model entity they discretize. Tags for elements and
        // nodes are globally unique (and not only per dimension, like entities).

        // A model entity of dimension 0 (a geometrical point) will contain a mesh
        // element of type point, as well as a mesh node. A model curve will contain
        // line elements as well as its interior nodes, while its boundary nodes
        // will be stored in the bounding model points. A model surface will contain
        // triangular and/or quadrangular elements and all the nodes not classified
        // on its boundary or on its embedded entities. A model volume will contain
        // tetrahedra, hexahedra, etc. and all the nodes not classified on its
        // boundary or on its embedded entities.

        // Get the mesh nodes for the entity (dim, tag):
        std::vector<std::size_t> nodeTags;
        std::vector<double> nodeCoords, nodeParams;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

        unsigned vert_offset = vertices.size();
        vertices.resize(vert_offset + nodeTags.size());
        for (unsigned i = 0; i < nodeTags.size(); i++)
        {
            vertices[vert_offset + i][0] = static_cast<Real>(nodeCoords[i*3]);
            vertices[vert_offset + i][1] = static_cast<Real>(nodeCoords[i*3 + 1]);
            vertices[vert_offset + i][2] = static_cast<Real>(nodeCoords[i*3 + 2]);
        }

        // Get the mesh elements for the entity (dim, tag):
        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

        // Extract all tetrahedra into a flat vector
        std::vector<unsigned> tetrahedra_vertex_indices;
        // Extract all surface triangles into a flat vector
        std::vector<unsigned> triangle_vertex_indices;
        for (unsigned i = 0; i < elemTypes.size(); i++)
        {
            if (elemTypes[i] == 4)
            {
                tetrahedra_vertex_indices.insert(tetrahedra_vertex_indices.end(), elemNodeTags[i].begin(), elemNodeTags[i].end());
            }

            if (elemTypes[i] == 2)
            {
                triangle_vertex_indices.insert(triangle_vertex_indices.end(), elemNodeTags[i].begin(), elemNodeTags[i].end());
            }
        }

        unsigned elem_offset = elements.size();
        unsigned num_tetrahedra = tetrahedra_vertex_indices.size()/4;
        elements.resize(elem_offset + num_tetrahedra);
        for (unsigned i = 0; i < num_tetrahedra; i++)
        {
            elements[elem_offset + i][0] = tetrahedra_vertex_indices[i*4] - 1;
            elements[elem_offset + i][1] = tetrahedra_vertex_indices[i*4 + 1] - 1;
            elements[elem_offset + i][2] = tetrahedra_vertex_indices[i*4 + 2] - 1;
            elements[elem_offset + i][3] = tetrahedra_vertex_indices[i*4 + 3] - 1;
        }

        unsigned face_offset = faces.size();
        unsigned num_faces = triangle_vertex_indices.size()/3;
        faces.resize(face_offset + num_faces);
        for (unsigned i = 0; i < num_faces; i++)
        {
            faces[face_offset + i][0] = triangle_vertex_indices[i*3] - 1;
            faces[face_offset + i][1] = triangle_vertex_indices[i*3 + 1] - 1;
            faces[face_offset + i][2] = triangle_vertex_indices[i*3 + 2] - 1;
        }
    }

    ParticleTetMesh tet_mesh(vertices, faces, elements);

    // write the loaded surface mesh part to file
    const std::string surface_mesh_filename = filename.substr(0,filename.length()-4) + "_surface_mesh.obj";
    tet_mesh.writeMeshToObjFile(surface_mesh_filename);

    return tet_mesh;

}

void MeshUtils::convertToSTL(const std::string& filename)
{
    std::cout << "MeshUtils::convertToSTL - converting " << filename << " to .stl format..." << std::endl;
    std::filesystem::path file_path(filename);

    Assimp::Importer importer;

    // And have it read the given file with some example postprocessing
    // Usually - if speed is not the most important aspect for you - you'll
    // probably to request more postprocessing than we do in this example.
    const aiScene* scene = importer.ReadFile( filename,
        aiProcess_Triangulate            |
        aiProcess_JoinIdenticalVertices  |
        aiProcess_FixInfacingNormals |
        aiProcess_SortByPType);

    // If the import failed, report it
    if (scene == nullptr)
    {
        std::cerr << "\tAssimp::Importer could not open " << filename << std::endl;
        std::cerr << "\tEnsure that the file is in a format that assimp can handle." << std::endl;
        return;
    }

    std::cout << "\tAssimp import successful" << std::endl;

    Assimp::Exporter exporter;
    const std::string& stl_filename = file_path.replace_extension(".stl").string();
    exporter.Export(scene, "stl", stl_filename);

    std::cout << "\tAssimp export successful - new file is " << stl_filename << "\n" << std::endl;
}

void MeshUtils::convertSTLtoMSH(const std::string& filename)
{
    std::cout << "MeshUtils::convertSTLtoMSH - converting " << filename << " from .stl to .msh format..." << std::endl;

    // ensure the file exists
    if (!std::filesystem::exists(filename))
    {
        std::cerr << "\t" << filename << " does not exist!" << std::endl;
        return;
    }


    std::filesystem::path file_path(filename);

    // ensure the file is an STL file
    if (file_path.extension() != ".stl" && file_path.extension() != ".STL")
    {
        std::cerr << "\t" << filename << " is not a .stl file!" << std::endl;
        return;
    }

    gmsh::open(filename);
    gmsh::option::setNumber("General.Verbosity", 5);

    // Get the surface (should be tag 1)
    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
    std::cout << "Surface tag: " << surfaces[0].second << std::endl;

    // Method 1: Force create volume with the actual surface tag
    int surface_tag = surfaces[0].second;
    int surface_loop_tag = gmsh::model::geo::addSurfaceLoop(std::vector<int>{surface_tag});
    int volume_tag = gmsh::model::geo::addVolume(std::vector<int>{surface_loop_tag});
    gmsh::model::geo::synchronize();

    std::cout << "Created volume with tag: " << volume_tag << std::endl;

    // Generate mesh
    gmsh::model::mesh::generate(3);


    const std::string& msh_filename = file_path.replace_extension(".msh").string();
    gmsh::write(msh_filename);

    std::cout << "\tGMSH conversion to .msh successful - new file is " << msh_filename << "\n" << std::endl;
    
}