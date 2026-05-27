#include "simobject/TetMeshObject.hpp"

#include "common/mesh/MeshUtils.hpp"

namespace SimObject
{

TetMeshObject::TetMeshObject(const Config::TetMeshObjectConfig& config)
    : Object_Base(config),
    _filename(config.filename())
{

}

void TetMeshObject::setup()
{
    // load mesh from file
    _mesh = MeshUtils::loadTetMeshFromGmshFile(_filename);

    // hard-coded material for now
    Real _E = 1e6;
    Real _nu = 0.45;
    Real _mu = _E / (2 * (1 + _nu));
    Real _lambda = (_E*_nu) / ( (1 + _nu) * (1 - 2*_nu) );

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
        }
    }
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

} // namespace SimObject