#include "simobject/TetMeshObject.hpp"

#include "common/mesh/MeshUtils.hpp"
#include "simulation/SimulationContext.hpp"

namespace SimObject
{

TetMeshObject::TetMeshObject(Sim::SimulationContext* ctx, const Config::TetMeshObjectConfig& config)
    : Object_Base(ctx, config),
    _filename(config.filename()),
    _E(config.E()),
    _nu(config.nu()),
    _density(config.density()),
    _config(config)
{

}

void TetMeshObject::setup()
{
    // load mesh from file
    _mesh = MeshUtils::loadTetMeshFromGmshFile(_filename, _ctx->particles);
    
    // scale the mesh according to the requested scaling
    _mesh.scale(_config.scaling());

    // rotate the mesh to the initial rotation
    Mat3r rot_mat = Math::RotMatFromXYZEulerAngles(_config.initialRotation());
    _mesh.rotateAbout(Vec3r::Zero(), rot_mat);

    // move the mesh to the initial position
    _mesh.moveTogether(_config.initialPosition());

    // set the current mesh state as the undeformed state
    _mesh.setCurrentStateAsUndeformedState();

    Real mu = _E / (2 * (1 + _nu));
    Real lambda = (_E*_nu) / ( (1 + _nu) * (1 - 2*_nu) );

    for (const auto& v_idx : _mesh.vertices())
    {
        // initialize masses to 0
        _ctx->particles.masses[v_idx] = 0;
    }

    // create Neo-Hookean energies for each element
    for (const auto& elem : _mesh.elements())
    {
        // get indices of vertices in the particle pool
        Vec4u pool_indices;
        for (int i = 0; i < 4; i++)
        {
            pool_indices[i] = _mesh.vertices().at(elem[i]);
        }

        // compute inverse undeformed basis
        const Vec3r& v1 = _mesh.vertex(elem[0]);
        const Vec3r& v2 = _mesh.vertex(elem[1]);
        const Vec3r& v3 = _mesh.vertex(elem[2]);
        const Vec3r& v4 = _mesh.vertex(elem[3]);
        Mat3r X;
        X.col(0) = (v1 - v4);
        X.col(1) = (v2 - v4);
        X.col(2) = (v3 - v4);
        Mat3r Q = X.inverse();
        Real rest_volume = std::abs(X.determinant() / 6.0);

        // create the new energy
        /** TODO: should I store this? */
        _ctx->energies.neo_hookean.addEnergy(
            pool_indices,
            lambda,
            mu,
            _config.dampingCoefficient(),
            Q,
            rest_volume
        );

        // add mass to each particle in this element
        for (int i = 0; i < 4; i++)
        {
            _ctx->particles.masses[pool_indices[i]] += 0.25 * rest_volume * _density;
        }
    }

    // add ground collision constraints for each particle
    for (auto& v_idx : _mesh.vertices())
    {
        _ctx->energies.ground_collision.addEnergy(v_idx);
    }


    // int fix_vertex = _mesh.face(0)[0];
    // auto& new_constraint = _attachment_constraints.emplace_back(&_mesh.particle(fix_vertex), _mesh.particle(fix_vertex).position);
    // auto& new_energy = _constraint_energies.emplace_back(&new_constraint, 1e1);
    // _mesh.particle(fix_vertex).energies.emplace_back(std::make_pair(&new_energy, 0));

    
    // Real total_mass = 0;
    // for (const auto& vertex : _mesh.vertices())
    // {
    //     total_mass += vertex.mass;
    // }
    // std::cout << "Total mass: " << total_mass << std::endl;
}

// void TetMeshObject::for_each_particle(std::function<void(Particle*)> func)
// {
//     for (auto& particle : _mesh.vertices())
//     {
//         func(&particle);
//     }
// }

// void TetMeshObject::for_each_particle(std::function<void(const Particle*)> func) const
// {
//     for (const auto& particle : _mesh.vertices())
//     {
//         func(&particle);
//     }
// }

// void TetMeshObject::for_each_energy(std::function<void(Energy::Energy_Base*)> func)
// {
//     for (auto& energy : _element_energies)
//         func(&energy);

//     for (auto& energy : _collision_energies)
//         func(&energy);

//     for (auto& energy : _constraint_energies)
//         func(&energy);
// }

// void TetMeshObject::for_each_energy(std::function<void(const Energy::Energy_Base*)> func) const
// {
//     for (const auto& energy : _element_energies)
//         func(&energy);
        
//     for (const auto& energy : _collision_energies)
//         func(&energy);

//     for (const auto& energy : _constraint_energies)
//         func(&energy);
// }

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