#include "simobject/Rod.hpp"

#include "simulation/SimulationContext.hpp"

#include "common/Math.hpp"

namespace SimObject
{

Rod::Rod(Sim::SimulationContext* ctx, const Config::RodConfig& config)
    : Object_Base(ctx, config)
    , _num_elements(config.elements())
    , _length(config.length())
    , _curvature(config.curvature())
    , _E(config.E())
    , _nu(config.nu())
    , _density(config.density())
    , _radius(config.diameter()/2)
{
    // compute shear modulus
    _G = _E / (2 * (1+_nu));

    // create particles
    int num_nodes = config.elements() + 1;
    _nodes.resize(num_nodes);

    Vec3r cur_position = config.initialPosition();
    Quaternion cur_rotation = Math::QuaternionFromXYZEulerAngles(config.initialRotation()); 

    Real ds = _length / (num_nodes - 1);
    Vec3r dR = _curvature * ds;
    for (int i = 0; i < num_nodes; i++)
    {
        _nodes[i] = _ctx->particles.addOrientedParticle(cur_position, cur_rotation, 0, Vec3r::Zero());
        _ctx->particles.velocities[_nodes[i]] = config.initialVelocity();

        
        cur_position += cur_rotation * Vec3r(0, 0, ds);
        cur_rotation = cur_rotation * Math::Exp_s3(dR);
    }

    // compute cross section properties
    _area = M_PI * _radius * _radius;
    _Ix = M_PI * _radius * _radius * _radius * _radius / 4.0;
    _Iz = 2*_Ix;

    // assign masses
    Real element_rest_length = _length / (num_nodes - 1);
    Real total_element_mass = element_rest_length * _area * _density;
    Vec3r total_element_rot_inertia = _density * element_rest_length * Vec3r(_Ix, _Iy, _Iz);
    for (unsigned e = 0; e < _num_elements; e++)
    {
        _ctx->particles.masses[_nodes[e]] += 0.5 * total_element_mass;
        _ctx->particles.rotationalInertia(_nodes[e+1]) += 0.5 * total_element_rot_inertia;
    }

    // create Cosserat energies
    for (unsigned e = 0; e < _num_elements; e++)
    {
        _ctx->energies.cosserat_rod.addEnergy(
            {_nodes[e], _nodes[e+1]},
            Vec6r(_G*_area, _G*_area, _E*_area, _E*_Ix, _E*_Iy, _G*_Iz),
            ds,
            _curvature
        );
    }

    // add fixed constraints
    if (config.baseFixed())
    {
        _ctx->particles.fixed[_nodes.front()] = true;
    }
    // if (config.baseFixed())
    // {
    //     _ctx->energies.one_sided_fixed_joint.addEnergy(
    //         _nodes.front(),
    //         Vec3r::Zero(),
    //         Quaternion::Identity(),
    //         _ctx->particles.positions[_nodes.front()],
    //         _ctx->particles.rotation(_nodes.front())
    //     );

    //     std::cout << "Adding fixed base constraint!" << std::endl;
    // }

    // if (config.tipFixed())
    // {
    //     _ctx->energies.one_sided_fixed_joint.addEnergy(
    //         _nodes.back(),
    //         Vec3r::Zero(),
    //         Quaternion::Identity(),
    //         _ctx->particles.positions[_nodes.back()],
    //         _ctx->particles.rotation(_nodes.back())
    //     );
    // }
    
}

void Rod::setup()
{

}

const Vec3r& Rod::nodePosition(unsigned idx) const
{
    return _ctx->particles.positions[_nodes[idx]];
}

const Quaternion& Rod::nodeRotation(unsigned idx) const
{
    return _ctx->particles.rotation(_nodes[idx]);
}

} // namespace SimObject