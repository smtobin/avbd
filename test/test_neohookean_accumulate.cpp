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
    Sim::SimulationContext ctx(1000, 0, 5000);

    // load ParticleTetMesh
    ParticleTetMesh mesh = MeshUtils::loadTetMeshFromGmshFile("../resource/cube2.msh", ctx.particles);

    // create NeoHookeanEnergy
    const Vec4u& indices = mesh.element(0);
    const Vec3r& v1 = mesh.vertex(indices[0]);
    const Vec3r& v2 = mesh.vertex(indices[1]);
    const Vec3r& v3 = mesh.vertex(indices[2]);
    const Vec3r& v4 = mesh.vertex(indices[3]);
    std::cout << "v1: " << v1.transpose() << std::endl;
    std::cout << "v2: " << v2.transpose() << std::endl;
    std::cout << "v3: " << v3.transpose() << std::endl;
    std::cout << "v4: " << v4.transpose() << std::endl;
    Mat3r X;
    X.col(0) = v1 - v4;
    X.col(1) = v2 - v4;
    X.col(2) = v3 - v4;
    Mat3r Q = X.inverse();

    unsigned e_idx = ctx.energies.neo_hookean.addEnergy(indices, 1e7, 1e6, 10, Q, 5);
    Energy::NeoHookeanEnergySolver::updateAfterTimeStep(
        e_idx,
        ctx.energies.neo_hookean,
        ctx.particles
    );

    // move the particles around a bit
    mesh.setVertex(indices[0], v1+Vec3r(0.3,0.3,0.3));
    mesh.setVertex(indices[1], v2+Vec3r(0.1,0.2,0.3));
    mesh.setVertex(indices[2], v3+Vec3r(-0.2,-0.3,-0.2));
    mesh.setVertex(indices[3], v4+Vec3r(-0.1,0.1,0.1));

    Real dt = 1e-3;
    for (unsigned l_idx = 0; l_idx < 4; l_idx++)
    {
        Mat3r H_single = Mat3r::Zero();
        Vec3r G_single = Vec3r::Zero();
        Energy::NeoHookeanEnergySolver::accumulate(
            e_idx,
            ctx.energies.neo_hookean,
            ctx.particles,
            l_idx,
            H_single,
            G_single,
            dt
        );

        Mat3r H4 = Mat3r::Zero();
        Vec3r G4 = Vec3r::Zero();
        unsigned e_idx4[4] = {e_idx, e_idx, e_idx, e_idx};
        unsigned l_idx4[4] = {l_idx, l_idx, l_idx, l_idx};
        Energy::NeoHookeanEnergySolver::accumulate4(
            e_idx4,
            ctx.energies.neo_hookean,
            ctx.particles,
            l_idx4,
            H4,
            G4,
            dt
        );

        std::cout << "G single: " << G_single.transpose() << std::endl;
        std::cout << "G AVX: " << G4.transpose() / 4 << std::endl;

        std::cout << "H single:\n" << H_single << std::endl;
        std::cout << "H AVX:\n" << H4/4 << std::endl;
    }

    /** Use finite differences to check linearizations */
    std::cout << "\n === FINITE DIFFERENCE CHECK === " << std::endl;
    for (unsigned l_idx = 0; l_idx < 4; l_idx++)
    {
        Mat3r H_orig = Mat3r::Zero();
        Vec3r G_orig = Vec3r::Zero();
        Energy::NeoHookeanEnergySolver::accumulate(
            e_idx,
            ctx.energies.neo_hookean,
            ctx.particles,
            l_idx,
            H_orig,
            G_orig,
            dt
        );

        Real E_orig = Energy::NeoHookeanEnergySolver::energy(
            e_idx,
            ctx.energies.neo_hookean,
            ctx.particles,
            dt
        );
        std::cout << "E_orig: " << E_orig << std::endl;

        Vec3r G_fd = Vec3r::Zero();
        Mat3r H_fd = Mat3r::Zero();
        for (unsigned k = 0; k < 3; k++)
        {
            Real delta = 1e-6;
            Vec3r delta_vec = Vec3r::Zero();
            delta_vec[k] = delta;

            Vec3r vi_cur = mesh.vertex(indices[l_idx]);
            mesh.setVertex(indices[l_idx], vi_cur + delta_vec);

            Mat3r H_new = Mat3r::Zero();
            Vec3r G_new = Vec3r::Zero();
            Energy::NeoHookeanEnergySolver::accumulate(
                e_idx,
                ctx.energies.neo_hookean,
                ctx.particles,
                l_idx,
                H_new,
                G_new,
                dt
            );
            Real E_new = Energy::NeoHookeanEnergySolver::energy(
                e_idx,
                ctx.energies.neo_hookean,
                ctx.particles,
                dt
            );
            std::cout << "E_new: " << E_new << std::endl;

            mesh.setVertex(indices[l_idx], vi_cur);
            
            G_fd[k] = (E_new - E_orig) / delta;
            H_fd.col(k) = (G_new - G_orig) / delta;
        }

        Mat3r H_diff = H_fd - H_orig;
        std::cout << "Local idx " << l_idx << std::endl;
        std::cout << "G_orig: " << G_orig.transpose() << std::endl;
        std::cout << "G_fd: " << G_fd.transpose() << std::endl;
        std::cout << "H_orig:\n" << H_orig << std::endl;
        std::cout << "H_fd:\n" << H_fd << std::endl;
        std::cout << "Local idx " << l_idx << " H_diff:\n" << H_diff << std::endl;

    }

    

    

}