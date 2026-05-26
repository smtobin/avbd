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

template <typename EnergyType>
bool testElementEnergy(ParticleTetMesh* mesh, int element_index)
{
    EnergyType energy(mesh, element_index, 1000, 1000);

    std::vector<Vec3r> numerical_gradients = numericalEnergyGradients(&energy);

    for (int pi = 0; pi < energy.numParticles(); pi++)
    {
        Vec3r gradient = energy.gradient(pi);
        Vec3r diff = gradient - numerical_gradients[pi];

        if (diff.norm() > 1e-4)
        {
            std::cout << "Gradient mismatch for particle " << pi << ":\n  " << "Numerical gradient: " << numerical_gradients[pi].transpose() << "\n  " <<
                "Analytical gradient: " << gradient.transpose() << std::endl;
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

    testElementEnergy<Energy::TetElementEnergy>(&mesh, element_index);
}