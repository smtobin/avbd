#include "energy/CosseratRodEnergyPool.hpp"
#include "energy/CosseratRodEnergySolver.hpp"

#include "simulation/SimulationContext.hpp"


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

int main()
{

    Real dt = 1e-3;
    // create SimulationContext
    Sim::SimulationContext ctx(1000, 1000, 5000, 1000, 1000);

    // add two oriented particles for rod element
    unsigned n1 = ctx.particles.addOrientedParticle(Vec3r::Zero(), Quaternion::Identity(), 10, Vec3r::Ones());
    unsigned n2 = ctx.particles.addOrientedParticle(Vec3r(0.3, 0.4, 0.5), Math::QuaternionFromXYZEulerAngles(Vec3r(60,-30,20)), 10, Vec3r::Ones());

    unsigned e_idx = ctx.energies.cosserat_rod.addEnergy(
        { n1, n2 },
        Vec6r(2, 2, 3.5, 0.5, 0.5, 0.8),
        0.2,
        Vec3r::Zero()
    );

    /** Use finite differences to check linearizations */
    std::cout << "\n === FINITE DIFFERENCE CHECK === " << std::endl;

    
    Vec2u indices = ctx.energies.cosserat_rod.data[e_idx].particle_indices;
    for (unsigned l_idx = 0; l_idx < 2; l_idx++)
    {
        auto F = [&](const Vec6r& dx) {
            Vec3r& p_cur = ctx.particles.positions[indices[l_idx]];
            Vec3r p_orig = p_cur;
            Quaternion& q_cur = ctx.particles.rotation(indices[l_idx]);
            Quaternion q_orig = q_cur;

            p_cur += dx.head<3>();
            q_cur = q_cur * Math::Exp_s3(dx.tail<3>());

            Real E = Energy::CosseratRodEnergySolver::energy(e_idx, ctx.energies.cosserat_rod, ctx.particles, dt);

            p_cur = p_orig;
            q_cur = q_orig;

            return E;
        };
        Vec6r G_fd = Gradient_FD(F);
        Mat6r H_fd = Hessian_FD(F);
        
        Vec6r G_orig = Vec6r::Zero();
        Mat6r H_orig = Mat6r::Zero();
        Energy::CosseratRodEnergySolver::accumulate(
            e_idx,
            ctx.energies.cosserat_rod,
            ctx.particles,
            l_idx,
            H_orig,
            G_orig,
            dt
        );

        Mat6r H_diff = H_fd - H_orig;
        std::cout << "Local idx " << l_idx << std::endl;
        std::cout << "G_orig: " << G_orig.transpose() << std::endl;
        std::cout << "G_fd: " << G_fd.transpose() << std::endl;
        std::cout << "H_orig:\n" << H_orig << std::endl;
        std::cout << "H_fd:\n" << H_fd << std::endl;
        std::cout << "Local idx " << l_idx << " H_diff:\n" << H_diff << std::endl;

    }

    

    

}