#include "energy/TetElementEnergy.hpp"

namespace Energy
{

TetElementEnergy::TetElementEnergy(const ParticleTetMesh* mesh, int element_index, Real lambda, Real mu, Real kd)
    : _mesh(mesh), _element(element_index), _lambda(lambda), _mu(mu), _kd(kd)
{

}

const Particle* TetElementEnergy::particle(int index) const
{
    const Vec4i& elem = _mesh->element(_element);
    return &_mesh->particle(elem[index]);
}

Real TetElementEnergy::energy(Real dt) const
{
    Real gamma = 1 + _mu / _lambda; 

    Mat3r F = _mesh->elementDeformationGradient(_element);
    Mat3r F_prev = _mesh->elementPreviousDeformationGradient(_element);

    Real V = _mesh->elementRestVolume(_element);

    Real detF_min1 = F.determinant() - gamma;
    Mat3r FTF = F.transpose() * F;
    Real strain_energy =  V * (_lambda/2 * detF_min1 * detF_min1 + _mu/2 * (FTF.trace() - 3));

    Mat3r E = FTF - Mat3r::Identity();
    Mat3r E_prev = F_prev.transpose() * F_prev - Mat3r::Identity();
    Mat3r E_dot = 1/dt * (E - E_prev);
    Real damping_energy = 0.5 * V * _kd * (E_dot.transpose() * E_dot).trace();

    return strain_energy + damping_energy;
}

Vec3r TetElementEnergy::gradient(int index, Real dt) const
{
    Mat3r F = _mesh->elementDeformationGradient(_element);
    Mat3r F_prev = _mesh->elementPreviousDeformationGradient(_element);
    Mat3r Q = _mesh->elementInvUndeformedBasis(_element);
    Real V = _mesh->elementRestVolume(_element);

    // gradient of hydrostatic part
    Real gamma = 1 + _mu / _lambda; 

    Mat3r F_cross;
    F_cross.col(0) = F.col(1).cross(F.col(2));
    F_cross.col(1) = F.col(2).cross(F.col(0));
    F_cross.col(2) = F.col(0).cross(F.col(1));

    Mat3r hyd_grad_full = V*_lambda * (F.determinant() - gamma) * F_cross * Q.transpose();
    
    // gradient of deviatoric part
    Mat3r dev_grad_full = V*_mu * F * Q.transpose();

    // gradient of damping energy
    Mat3r E = F.transpose() * F - Mat3r::Identity();
    Mat3r E_prev = F_prev.transpose() * F_prev - Mat3r::Identity();
    Mat3r E_dot = 1/dt * (E - E_prev);

    if (index < 3)
    {
        return hyd_grad_full.col(index) + dev_grad_full.col(index) + 2/dt * V * _kd * F * E_dot * Q.row(index).transpose();
    }
    else
    {
        Vec3r hyd_grad = -hyd_grad_full.col(0) - hyd_grad_full.col(1) - hyd_grad_full.col(2);
        Vec3r dev_grad = -dev_grad_full.col(0) - dev_grad_full.col(1) - dev_grad_full.col(2);
        return hyd_grad + dev_grad - 2/dt * V * _kd * F * E_dot * ( Q.row(0).transpose() + Q.row(1).transpose() + Q.row(2).transpose() );
    }
    
}

Mat3r TetElementEnergy::hessian(int index, Real dt) const
{
    Mat3r F = _mesh->elementDeformationGradient(_element);
    Mat3r F_prev = _mesh->elementPreviousDeformationGradient(_element);
    Mat3r Q = _mesh->elementInvUndeformedBasis(_element);
    Real V = _mesh->elementRestVolume(_element);

    // hessian of hydrostatic part
    Mat3r F_cross;
    F_cross.col(0) = F.col(1).cross(F.col(2));
    F_cross.col(1) = F.col(2).cross(F.col(0));
    F_cross.col(2) = F.col(0).cross(F.col(1));

    Mat3r detF_grad_full = F_cross * Q.transpose();

    Mat3r E = F.transpose() * F - Mat3r::Identity();
    Mat3r E_prev = F_prev.transpose() * F_prev - Mat3r::Identity();
    Mat3r E_dot = 1/dt * (E - E_prev);
    
    if (index < 3)
    {
        Mat3r hyd_hess = V*_lambda * detF_grad_full.col(index) * detF_grad_full.col(index).transpose();
        Mat3r dev_hess = V*_mu * Q.row(index).squaredNorm() * Mat3r::Identity();

        Vec3r qi = Q.row(index);
        Vec3r Fqi = F * qi;
        
        Mat3r damp_hess = 2/dt * V * _kd * ( (qi.transpose() * E_dot * qi) * Mat3r::Identity() + 1/dt * (Fqi * Fqi.transpose() + F*F.transpose() * qi.squaredNorm()) );

        return hyd_hess + dev_hess + damp_hess;
    }
    else
    {
        Vec3r a3 =
            -detF_grad_full.col(0)
            -detF_grad_full.col(1)
            -detF_grad_full.col(2);

        Vec3r q3 =
            -Q.row(0).transpose()
            -Q.row(1).transpose()
            -Q.row(2).transpose();

        Mat3r hyd_hess =
            V*_lambda * a3 * a3.transpose();

        Mat3r dev_hess =
            V*_mu * q3.squaredNorm() * Mat3r::Identity();

        Vec3r qi = Q.row(0) + Q.row(1) + Q.row(2);
        Vec3r Fqi = F * qi;
        
        Mat3r damp_hess = 2/dt * V * _kd * ( (qi.transpose() * E_dot * qi) * Mat3r::Identity() + 1/dt * (Fqi * Fqi.transpose() + F*F.transpose() * qi.squaredNorm()) );


        return hyd_hess + dev_hess + damp_hess;
    }
}

} // namespace Energy