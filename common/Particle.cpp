#include "common/Particle.hpp"
#include "energy/EnergyBase.hpp"

void Particle::inertialUpdate(Real dt, const Vec3r& a_ext)
{
    // use adaptive initialization (Sec 3.7 in VBD paper)
    Vec3r a = (velocity - prev_velocity) / dt;

    Real a_ext_norm = a_ext.norm();
    Real a_along_a_ext = a.dot(a_ext) / a_ext_norm;

    Real a_tilde;
    if (a_along_a_ext > a_ext_norm)
        a_tilde = 1;
    else if (a_along_a_ext < 0)
        a_tilde = 0;
    else
        a_tilde = a_along_a_ext / a_ext_norm;

    Vec3r a_tilde_vec = a_tilde * a_ext;

    // compute the inertially predicted position
    inertial_position = position + dt*velocity + dt*dt*a_ext;

    // move particle to its initialized position
    position += dt*velocity + dt*dt*a_tilde_vec;
}

void Particle::solveParticle(Real dt)
{
    Real kd = 1e-4;

    // iterate through energies, and sum up force and Hessian contributions
    Vec3r energy_force = Vec3r::Zero();
    Mat3r energy_hess = Mat3r::Zero();
    for (const auto& pair : energies)
    {
        const auto& energy = pair.first;
        int index = pair.second;


        Mat3r this_hess = energy->hessian(index);
        energy_force += energy->gradient(index) + kd/dt * this_hess * (position - prev_position);
        energy_hess += (1 + kd/dt) * this_hess;
    }

    // assemble LHS and RHS of single-particle system
    Vec3r RHS = -mass / (dt*dt) * (position - inertial_position) - energy_force;
    Mat3r LHS = mass / (dt*dt) * Mat3r::Identity() + energy_hess;

    Vec3r dx = LHS.partialPivLu().solve(RHS);

    position += dx;
}

void Particle::velocityUpdate(Real dt)
{
    prev_velocity = velocity;
    velocity = (position - prev_position) / dt;

    prev_position = position;
}

