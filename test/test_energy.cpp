#include "energy/EnergyBase.hpp"
#include "energy/TetElementEnergy.hpp"

#include "common/mesh/MeshUtils.hpp"

#include <gmsh.h>

std::vector<Vec3r> numericalEnergyGradients(const Energy::Energy_Base* energy)
{
    Real orig_energy = energy->energy();
    std::vector<Vec3r> gradients;
    
    Real delta = 1e-8;
    for (int pi = 0; pi < energy->numParticles(); pi++)
    {
        Particle* particle_i = const_cast<Particle*>(energy->particle(pi));
        Vec3r gradient = Vec3r::Zero();
        for (int k = 0; k < 3; k++)
        {
            particle_i->position[k] += delta;
            Real new_energy = energy->energy();
            particle_i->position[k] -= delta;

            gradient[k] = (new_energy - orig_energy) / delta;
        }

        gradients.push_back(gradient);
    }

    return gradients;
}

std::vector<Mat3r> numericalEnergyHessians(const Energy::Energy_Base* energy)
{
    Real delta = 1e-8;
    std::vector<Mat3r> hessians; 
    for (int pi = 0; pi < energy->numParticles(); pi++)
    {
        Vec3r orig_grad = energy->gradient(pi);

        Particle* particle_i = const_cast<Particle*>(energy->particle(pi));
        Mat3r hessian = Mat3r::Zero();
        for (int k = 0; k < 3; k++)
        {
            particle_i->position[k] += delta;
            Vec3r new_grad = energy->gradient(pi);
            particle_i->position[k] -= delta;

            hessian.col(k) = (new_grad - orig_grad) / delta;
        }

        hessians.push_back(hessian);
    }

    return hessians;
}

template <typename EnergyType>
bool testElementEnergy(ParticleTetMesh* mesh, int element_index)
{
    EnergyType energy(mesh, element_index, 1000, 1000);

    std::vector<Vec3r> numerical_gradients = numericalEnergyGradients(&energy);
    std::vector<Mat3r> numerical_hessians = numericalEnergyHessians(&energy);

    for (int pi = 0; pi < energy.numParticles(); pi++)
    {
        Vec3r gradient = energy.gradient(pi);
        Vec3r diff = gradient - numerical_gradients[pi];

        if (diff.norm() > 1e-4)
        {
            std::cout << "Gradient mismatch for particle " << pi << ":\n  " << "Numerical gradient: " << numerical_gradients[pi].transpose() << "\n  " <<
                "Analytical gradient: " << gradient.transpose() << std::endl;
        }

        Mat3r hessian = energy.hessian(pi);
        Mat3r diff_hess = hessian - numerical_hessians[pi];
        if (diff_hess.norm() > 1e-4)
        {
            std::cout << "Hessian mismatch for particle " << pi << ":\n  " << "Numerical hessian:\n" << numerical_hessians[pi] << "\n  " <<
                "Analytical hessian:\n" << hessian << std::endl;
        }
    }

    return false;
}

int main()
{
    gmsh::initialize();

    // load initial mesh
    ParticleTetMesh mesh = MeshUtils::loadTetMeshFromGmshFile("../resource/single.msh");
    int element_index = 0;

    // perturb the tet a little bit
    mesh.setVertex(0, Vec3r(-0.4,-0.2,-0.3));

    testElementEnergy<Energy::TetElementEnergy>(&mesh, element_index);
}