#pragma once

#include "common/mesh/ParticleMesh.hpp"
#include "common/mesh/ParticleTetMesh.hpp"

struct MeshUtils
{
    static ParticleMesh loadSurfaceMeshFromFile(const std::string& filename, ParticlePool& pool);

    static void convertToSTL(const std::string& filename);

    static void convertSTLtoMSH(const std::string& filename);

    static ParticleTetMesh loadTetMeshFromGmshFile(const std::string& filename, ParticlePool& pool);
};