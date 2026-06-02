#include "energy/TetElementEnergy.hpp"

namespace Energy
{

TetElementEnergy::TetElementEnergy(const ParticleTetMesh* mesh, int element_index, Real lambda, Real mu)
    : _mesh(mesh), _element(element_index), _lambda(lambda), _mu(mu)
{

}

const Particle* TetElementEnergy::particle(int index) const
{
    const Vec4i& elem = _mesh->element(_element);
    return &_mesh->particle(elem[index]);
}

Real TetElementEnergy::energy() const
{

    Real gamma = 1 + _mu / _lambda; 

    Mat3r F = _mesh->elementDeformationGradient(_element);
    Real V = _mesh->elementRestVolume(_element);

    Real detF_min1 = F.determinant() - gamma;
    Mat3r FTF = F.transpose() * F;
    return V * (_lambda/2 * detF_min1 * detF_min1 + _mu/2 * (FTF.trace() - 3));
}

Vec3r TetElementEnergy::gradient(int index) const
{
    Mat3r F = _mesh->elementDeformationGradient(_element);
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

    if (index < 3)
    {

        return hyd_grad_full.col(index) + dev_grad_full.col(index);
    }
    else
    {
        Vec3r hyd_grad = -hyd_grad_full.col(0) - hyd_grad_full.col(1) - hyd_grad_full.col(2);
        Vec3r dev_grad = -dev_grad_full.col(0) - dev_grad_full.col(1) - dev_grad_full.col(2);
        return hyd_grad + dev_grad;
    }
    
}

Mat3r TetElementEnergy::hessian(int index) const
{
    Mat3r F = _mesh->elementDeformationGradient(_element);
    Mat3r Q = _mesh->elementInvUndeformedBasis(_element);
    Real V = _mesh->elementRestVolume(_element);

    // hessian of hydrostatic part
    Mat3r F_cross;
    F_cross.col(0) = F.col(1).cross(F.col(2));
    F_cross.col(1) = F.col(2).cross(F.col(0));
    F_cross.col(2) = F.col(0).cross(F.col(1));

    Mat3r detF_grad_full = F_cross * Q.transpose();
    
    if (index < 3)
    {
        Mat3r hyd_hess = V*_lambda * detF_grad_full.col(index) * detF_grad_full.col(index).transpose();
        Mat3r dev_hess = V*_mu * Q.row(index).squaredNorm() * Mat3r::Identity();

        return hyd_hess + dev_hess;
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
            _lambda * a3 * a3.transpose();

        Mat3r dev_hess =
            _mu * q3.squaredNorm() * Mat3r::Identity();

        return hyd_hess + dev_hess;
    }
}

} // namespace Energy