#include "energy/joint/OneSidedFixedJointEnergySolver.hpp"

#include "simulation/SimulationContext.hpp"
#include "test/finite_difference.hpp"

int main()
{

    Real dt = 1e-3;
    // create SimulationContext
    Sim::SimulationContext ctx(1000, 1000, 5000, 1000, 1000);

    // add two oriented particles for rod element
    unsigned n1 = ctx.particles.addOrientedParticle(Vec3r::Zero(), Quaternion::Identity(), 10, Vec3r::Ones());

    unsigned e_idx = ctx.energies.one_sided_fixed_joint.addEnergy(
        n1,
        Vec3r(0.4, 0.6, 0.2),
        Math::QuaternionFromXYZEulerAngles(Vec3r(60,-30,20)),
        Vec3r(-0.3, -0.5, -0.2),
        Math::QuaternionFromXYZEulerAngles(Vec3r(-60,50,-20))
    );

    /** Use finite differences to check linearizations */
    std::cout << "\n === FINITE DIFFERENCE CHECK === " << std::endl;

    
    Vec1u indices = ctx.energies.one_sided_fixed_joint.data[e_idx].particle_indices;
    for (unsigned l_idx = 0; l_idx < 1; l_idx++)
    {
        auto F = [&](const Vec6r& dx) {
            Vec3r& p_cur = ctx.particles.positions[indices[l_idx]];
            Vec3r p_orig = p_cur;
            Quaternion& q_cur = ctx.particles.rotation(indices[l_idx]);
            Quaternion q_orig = q_cur;

            p_cur += dx.head<3>();
            q_cur = q_cur * Math::Exp_s3(dx.tail<3>());

            Real E = Energy::OneSidedFixedJointEnergySolver::energy(e_idx, ctx.energies.one_sided_fixed_joint, ctx.particles, dt);

            p_cur = p_orig;
            q_cur = q_orig;

            return E;
        };
        Vec6r G_fd = Gradient_FD(F);
        Mat6r H_fd = Hessian_FD(F);
        
        Vec6r G_orig = Vec6r::Zero();
        Mat6r H_orig = Mat6r::Zero();
        Energy::OneSidedFixedJointEnergySolver::accumulate(
            e_idx,
            ctx.energies.one_sided_fixed_joint,
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