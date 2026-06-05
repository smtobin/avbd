#include "energy/EnergyBase.hpp"
#include "energy/TetElementEnergy.hpp"
#include "energy/GroundCollisionEnergy.hpp"

#include "common/mesh/MeshUtils.hpp"

#include <gmsh.h>

#define DT 1e-3

std::vector<Vec3r> numericalEnergyGradients(const Energy::Energy_Base* energy)
{
    Real orig_energy = energy->energy(DT);
    std::vector<Vec3r> gradients;
    
    Real delta = 1e-8;
    for (int pi = 0; pi < energy->numParticles(); pi++)
    {
        Particle* particle_i = const_cast<Particle*>(energy->particle(pi));
        Vec3r gradient = Vec3r::Zero();
        for (int k = 0; k < 3; k++)
        {
            particle_i->position[k] += delta;
            Real new_energy = energy->energy(DT);
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
        Vec3r orig_grad = energy->gradient(pi, DT);

        Particle* particle_i = const_cast<Particle*>(energy->particle(pi));
        Mat3r hessian = Mat3r::Zero();
        for (int k = 0; k < 3; k++)
        {
            particle_i->position[k] += delta;
            Vec3r new_grad = energy->gradient(pi, DT);
            particle_i->position[k] -= delta;

            hessian.col(k) = (new_grad - orig_grad) / delta;
        }

        hessians.push_back(hessian);
    }

    return hessians;
}

bool testEnergy(const Energy::Energy_Base* energy)
{
    std::vector<Vec3r> numerical_gradients = numericalEnergyGradients(energy);
    std::vector<Mat3r> numerical_hessians = numericalEnergyHessians(energy);

    for (int pi = 0; pi < energy->numParticles(); pi++)
    {
        Vec3r gradient = energy->gradient(pi, DT);
        Vec3r diff = gradient - numerical_gradients[pi];

        if (diff.norm() > 1e-4)
        {
            std::cout << "Gradient mismatch for particle " << pi << ":\n  " << "Numerical gradient: " << numerical_gradients[pi].transpose() << "\n  " <<
                "Analytical gradient: " << gradient.transpose() << std::endl;
        }

        Mat3r hessian = energy->hessian(pi, DT);
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
    mesh.setCurrentStateAsUndeformedState();

    // std::cout << "F: " << mesh.elementDeformationGradient(0) << std::endl;

    // perturb the tet a little bit
    // std::cout << "Mesh vertex 0 pos: " << mesh.vertex(0).transpose() << std::endl;
    mesh.setVertex(0, Vec3r(-0.4,-0.2,-0.3));
    // std::cout << "Mesh vertex 0 pos: " << mesh.vertex(0).transpose() << std::endl;

    // std::cout << "F: " << mesh.elementDeformationGradient(0) << std::endl;

    Real _E = 1e5;
    Real _nu = 0.3;
    Real _mu = _E / (2 * (1 + _nu));
    Real _lambda = (_E*_nu) / ( (1 + _nu) * (1 - 2*_nu) );
    Real _kd = 1e-2;
    Energy::TetElementEnergy tet_energy(&mesh, element_index, _lambda, _mu, _kd);
    testEnergy(&tet_energy);

    // Energy::GroundCollisionEnergy ground_collision_energy(&mesh.particle(0));
    // testEnergy(&ground_collision_energy);


}