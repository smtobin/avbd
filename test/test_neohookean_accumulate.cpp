#include "energy/NeoHookeanEnergyPool.hpp"
#include "energy/NeoHookeanEnergySolver.hpp"

#include "simulation/SimulationContext.hpp"

#include "common/mesh/ParticleTetMesh.hpp"
#include "common/mesh/MeshUtils.hpp"

#include <gmsh.h>

int main()
{
    gmsh::initialize();
    
    // create SimulationContext
    Sim::SimulationContext ctx(1000, 5000);

    // load ParticleTetMesh
    ParticleTetMesh mesh = MeshUtils::loadTetMeshFromGmshFile("../resource/cube2.msh", ctx.particles);

    // create NeoHookeanEnergy
    const Vec4u& indices = mesh.element(0);
    const Vec3r& v1 = mesh.vertex(indices[0]);
    const Vec3r& v2 = mesh.vertex(indices[1]);
    const Vec3r& v3 = mesh.vertex(indices[2]);
    const Vec3r& v4 = mesh.vertex(indices[3]);
    Mat3r X;
    X.col(0) = v1 - v4;
    X.col(1) = v2 - v4;
    X.col(2) = v3 - v4;
    Mat3r Q = X.inverse();

    unsigned e_idx = ctx.energies.neo_hookean.addEnergy(indices, 1000, 1000, 100, Q, 5);
    Energy::NeoHookeanEnergySolver::updateAfterTimeStep(
        e_idx,
        ctx.energies.neo_hookean,
        ctx.particles
    );

    Mat3r H_single;
    Vec3r G_single;
    Energy::NeoHookeanEnergySolver::accumulate(
        e_idx,
        ctx.energies.neo_hookean,
        ctx.particles,
        3,
        H_single,
        G_single,
        1e-3
    );

    Mat3r H4;
    Vec3r G4;
    unsigned e_idx4[4] = {e_idx, e_idx, e_idx, e_idx};
    unsigned l_idx4[4] = {3, 3, 3, 3};
    Energy::NeoHookeanEnergySolver::accumulate4(
        e_idx4,
        ctx.energies.neo_hookean,
        ctx.particles,
        l_idx4,
        H4,
        G4,
        1e-3
    );

    std::cout << "G single: " << G_single.transpose() << std::endl;
    std::cout << "G AVX: " << G4.transpose() / 4 << std::endl;

}