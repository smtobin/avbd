#include "simobject/TetMeshObject.hpp"

#include "common/mesh/MeshUtils.hpp"

namespace SimObject
{

TetMeshObject::TetMeshObject(const Config::TetMeshObjectConfig& config)
    : Object_Base(config),
    _filename(config.filename()),
    _config(config)
{

}

void TetMeshObject::setup()
{
    // load mesh from file
    _mesh = MeshUtils::loadTetMeshFromGmshFile(_filename);
    _mesh.moveTo(_config.initialPosition());
    _mesh.setCurrentStateAsUndeformedState();

    // hard-coded material for now
    Real _E = 1e6;
    Real _nu = 0.3;
    Real _mu = _E / (2 * (1 + _nu));
    Real _lambda = (_E*_nu) / ( (1 + _nu) * (1 - 2*_nu) );

    Real _density = 1000;

    for (auto& vertex : _mesh.vertices())
    {
        vertex.mass = 0;
    }

    // create energies for each element
    _element_energies.reserve(_mesh.numElements());
    for (const auto& elem_index : _mesh.elements().validIndices())
    {
        auto& new_energy = _element_energies.emplace_back(&_mesh, elem_index, _lambda, _mu);

        // add newly created energy to each particle in this element
        const Vec4i& elem = _mesh.element(elem_index);
        for (int k = 0; k < 4; k++)
        {
            _mesh.particle(elem[k]).energies.emplace_back(std::make_pair(&new_energy, k));
            _mesh.particle(elem[k]).mass += 0.25*_mesh.elementVolume(elem_index) * _density;
        }
    }

    // add ground constraints for each particle
    _collision_energies.reserve(_mesh.numVertices());
    for (auto& vertex : _mesh.vertices())
    {
        auto& new_energy = _collision_energies.emplace_back(&vertex);
        vertex.energies.emplace_back(std::make_pair(&new_energy, 0));
    }

    
    // Real total_mass = 0;
    // for (const auto& vertex : _mesh.vertices())
    // {
    //     total_mass += vertex.mass;
    // }
    // std::cout << "Total mass: " << total_mass << std::endl;
}

void TetMeshObject::for_each_particle(std::function<void(Particle*)> func)
{
    for (auto& particle : _mesh.vertices())
    {
        func(&particle);
    }
}

void TetMeshObject::for_each_particle(std::function<void(const Particle*)> func) const
{
    for (const auto& particle : _mesh.vertices())
    {
        func(&particle);
    }
}

void TetMeshObject::for_each_energy(std::function<void(Energy::Energy_Base*)> func)
{
    for (auto& energy : _element_energies)
        func(&energy);

    for (auto& energy : _collision_energies)
        func(&energy);
}

void TetMeshObject::for_each_energy(std::function<void(const Energy::Energy_Base*)> func) const
{
    for (const auto& energy : _element_energies)
        func(&energy);
        
    for (const auto& energy : _collision_energies)
        func(&energy);
}

/** Provides a way to iterate through all QUADRATIC energies owned by the object. */
// void TetMeshObject::for_each_quadratic_energy(std::function<void(QuadraticEnergy*)> func)
// {       
//     for (auto& energy : _collision_energies)
//         func(&energy);
// }

// void TetMeshObject::for_each_quadratic_energy(std::function<void(QuadraticEnergy*)> func) const
// {
//     for (const auto& energy : _collision_energies)
//         func(&energy);
// }

} // namespace SimObject