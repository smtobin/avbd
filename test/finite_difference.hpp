#pragma once

#include "common/common.hpp"

/** Function F(dx) is a scalar function that takes combined perturbation [dp, dR] and outputs the new scalar output. */
template<typename Func>
Mat6r Hessian_FD(Func&& F)
{
    Mat6r H;
    Vec6r ei, ej;
    Real eps = 1e-6;
    for (unsigned i = 0; i < 6; i++)
    {
        ei = Vec6r::Zero();
        ei[i] = 1;
        for (unsigned j = 0; j < 6; j++)
        {
            ej = Vec6r::Zero();
            ej[j] = 1;

            if (i == j)
            {
                H(i,i) = ( F(eps*ei) - 2*F(Vec6r::Zero()) + F(-eps*ei) ) / (eps*eps);
            }
            else
            {
                H(i,j) = ( F(eps*(ei + ej)) - F(eps*(ei - ej)) - F(eps*(-ei + ej)) + F(-eps*(ei + ej)) ) / (4*eps*eps);
            }
        }
    }

    return H;
}

template<typename Func>
Vec6r Gradient_FD(Func&& F)
{
    Vec6r G;
    Vec6r ei;
    Real eps = 1e-6;
    for (unsigned i = 0; i < 6; i++)
    {
        ei = Vec6r::Zero();
        ei[i] = 1;

        G[i] = ( F(eps*ei) - F(Vec6r::Zero()) ) / eps;
    }

    return G;
}