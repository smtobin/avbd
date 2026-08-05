#pragma once

#include "common/common.hpp"
#include "common/AVX.hpp"
#include "common/ParticlePool.hpp"
#include "energy/NeoHookeanEnergyPool.hpp"

namespace Energy
{

/** Implements the stable Neo-Hookean per-element energies, seen in Macklin et al. 2021 */
struct NeoHookeanEnergySolver
{
    /** Public typedefs */
    using PoolType = NeoHookeanEnergyPool;
    static constexpr bool SupportsPositional = true;
    static constexpr bool SupportsOriented = false;

    /** Energy functional */
    static Real energy(
        unsigned e_idx,
        const NeoHookeanEnergyPool& energies,
        ParticlePool& particles,
        Real dt
    )
    {
        const Vec4u& indices = energies.data[e_idx].particle_indices;
        const Mat3r& Q = energies.data[e_idx].Q;
        const Real V = energies.data[e_idx].rest_volume;
        const Real kd = energies.data[e_idx].kd;
        const Real mu = energies.data[e_idx].mu;
        const Real lambda = energies.data[e_idx].lambda;

        Mat3r F;
        computeF(
            particles.positions[indices[0]],
            particles.positions[indices[1]],
            particles.positions[indices[2]],
            particles.positions[indices[3]],
            Q,
            F
        );

        Mat3r E = 0.5*( F.transpose()*F - Mat3r::Identity() );
        Mat3r Edot = 1/dt*(E - energies.data[e_idx].E_prev);
        Real damp_energy = V * kd * Edot.squaredNorm();

        Real detF_fac = (F.determinant() - mu/lambda - 1);
        Real hyd_energy = 0.5 * lambda * V * detF_fac * detF_fac;
        Real dev_energy = 0.5 * mu * V * ( (F.transpose()*F).trace() - 3);

        return hyd_energy + dev_energy + damp_energy;
    }

    /** Required - does nothing */
    static void updateAfterIteration(
        unsigned /* e_idx */,
        const NeoHookeanEnergyPool& /* energies */,
        ParticlePool& /* particles */
    )
    {

    }

    /** Computes the previous deformation gradient F and the previous Green strain E
     * @param e_idx : the energy index
     * @param energies : the memory pool for the energy
     * @param particles : the simulation particle memory pool
     */
    static void updateAfterTimeStep(
        unsigned e_idx,
        NeoHookeanEnergyPool& energies,
        ParticlePool& particles
    )
    {
        const Vec4u& indices = energies.data[e_idx].particle_indices;
        const Mat3r& Q = energies.data[e_idx].Q;

        Mat3r F_prev;
        computeF(
            particles.previous_positions[indices[0]],
            particles.previous_positions[indices[1]],
            particles.previous_positions[indices[2]],
            particles.previous_positions[indices[3]],
            Q,
            F_prev
        );

        energies.data[e_idx].E_prev = 0.5 * (F_prev.transpose() * F_prev - Mat3r::Identity());
    }

    /** Helper function to compute the deformation gradient F given 4 vertices and the rest-state matrix Q.
     * 
     * F = [v1 - v4   v2 - v4   v3 - v4] * Q
     */
    static void computeF(
        const Vec3r& v1,
        const Vec3r& v2,
        const Vec3r& v3, 
        const Vec3r& v4,
        const Mat3r& Q,
        Mat3r& F)
    {
        F.col(0).noalias() = v1 - v4;
        F.col(1).noalias() = v2 - v4;
        F.col(2).noalias() = v3 - v4;

        F = F*Q;
    }

    /** Computes the Hessian and gradient for a specified particle affected by this energy.
     * Updates the accumulated vertex Hessian and gradients.
     * @param e_idx : the energy index
     * @param energies : the memory pool for the energies
     * @param particles : the memory pool for the particles
     * @param local_idx : the "local" index of the particle in the energy
     * @param particle_H : the current accumulated particle Hessian - this is updated by this function
     * @param particle_G : the current accumulated particle gradient - this is updated by this function
     * @param dt : the time step
     */
    static void accumulate(
        unsigned e_idx,
        const NeoHookeanEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx,
        Mat3r& particle_H,
        Vec3r& particle_G,
        Real dt
    )
    {
        // extract data from energy pool
        // const Vec4u& indices = energies.particle_indices[e_idx];
        // Real V = energies.rest_volumes[e_idx];
        // const Mat3r& Q = energies.Qs[e_idx];
        // Real lambda = energies.lambdas[e_idx];
        // Real mu = energies.mus[e_idx];
        // Real kd = energies.kds[e_idx];

        const Vec4u& indices = energies.data[e_idx].particle_indices;
        Real V = energies.data[e_idx].rest_volume;
        const Mat3r& Q = energies.data[e_idx].Q;
        Real lambda = energies.data[e_idx].lambda;
        Real mu = energies.data[e_idx].mu;
        Real kd = energies.data[e_idx].kd;
        const Mat3r& E_prev = energies.data[e_idx].E_prev;

        // compute F for this timestep and the previous timestep
        Mat3r F;
        computeF(
            particles.positions[indices[0]],
            particles.positions[indices[1]],
            particles.positions[indices[2]],
            particles.positions[indices[3]],
            Q, 
            F
        );

        // hydrostatic gradient
        Vec3r c0 = F.col(1).cross(F.col(2));
        Vec3r c1 = F.col(2).cross(F.col(0));
        Vec3r c2 = F.col(0).cross(F.col(1));

        // Mat3r F_cross;
        // F_cross.col(0) = F.col(1).cross(F.col(2));
        // F_cross.col(1) = F.col(2).cross(F.col(0));
        // F_cross.col(2) = F.col(0).cross(F.col(1));

        // Mat3r detF_grad_full = F_cross * Q.transpose();
        // Mat3r hyd_grad_full = V*lambda * (F.determinant() - mu/lambda - 1) * detF_grad_full;

        // deviatoric gradient
        // Mat3r dev_grad_full = V*mu * F * Q.transpose();

        const Real detF = F.col(0).dot(c0);

        

        // compute the Hessian and gradient
        Vec3r qi;
        Real sign;
        if (local_idx < 3)
        {
            qi = Q.row(local_idx);
            sign = 1.0;
        }
        else
        {
            qi = Q.row(0) + Q.row(1) + Q.row(2);
            sign = -1.0;
        }

        const Vec3r Fqi = F*qi;
        const Vec3r Fcross_qi = c0*qi[0] + c1*qi[1] + c2*qi[2];
        
        const Real hyd_mult = V*lambda * (detF - mu/lambda - 1);
        const Real dev_mult = V*mu;

        // compute rate of change of Green strain
        Mat3r FtF;
        FtF.noalias() = F.transpose() * F;
        
        const Real qi_norm2 = qi.squaredNorm();

        // gradient
        Real damp_mult_factor = (2/dt * V * kd);
        particle_G += sign * hyd_mult * Fcross_qi;
        particle_G += sign * dev_mult * Fqi;

        Mat3r E_dot = 0.5 * FtF;
        E_dot.diagonal().array() -= 0.5;
        E_dot -= E_prev;
        E_dot *= 1/dt;

        particle_G += (sign * damp_mult_factor) * (F * (E_dot * qi));

        // Hessian
        // hydrostatic
        particle_H += (V * lambda) * Fcross_qi * Fcross_qi.transpose();
        // deviatoric
        particle_H.diagonal().array() += V * mu * qi_norm2;

        // damping
        Real qi_Edot_qi = qi.transpose() * E_dot * qi;
        particle_H.diagonal().array() += damp_mult_factor * qi_Edot_qi;

        particle_H += (0.5 * damp_mult_factor/dt) * (F*F.transpose() * qi_norm2 + Fqi * Fqi.transpose());
    }

    struct Batch4
    {
        alignas(32) Real qx[4];
        alignas(32) Real qy[4];
        alignas(32) Real qz[4];

        alignas(32) Real sign[4];

        alignas(32) Real V[4];
        alignas(32) Real lambda[4];
        alignas(32) Real mu[4];
        alignas(32) Real kd[4];

        alignas(32) Real p0x[4];
        alignas(32) Real p0y[4];
        alignas(32) Real p0z[4];

        alignas(32) Real p1x[4];
        alignas(32) Real p1y[4];
        alignas(32) Real p1z[4];

        alignas(32) Real p2x[4];
        alignas(32) Real p2y[4];
        alignas(32) Real p2z[4];

        alignas(32) Real p3x[4];
        alignas(32) Real p3y[4];
        alignas(32) Real p3z[4];

        alignas(32) Real q00[4];
        alignas(32) Real q01[4];
        alignas(32) Real q02[4];
        alignas(32) Real q10[4];
        alignas(32) Real q11[4];
        alignas(32) Real q12[4];
        alignas(32) Real q20[4];
        alignas(32) Real q21[4];
        alignas(32) Real q22[4];

        alignas(32) Real Eprev00[4];
        alignas(32) Real Eprev01[4];
        alignas(32) Real Eprev02[4];
        alignas(32) Real Eprev10[4];
        alignas(32) Real Eprev11[4];
        alignas(32) Real Eprev12[4];
        alignas(32) Real Eprev20[4];
        alignas(32) Real Eprev21[4];
        alignas(32) Real Eprev22[4];
    };

    static void computeF4(
        const Batch4& b,
        AVX::Mat3x3Packet& F
    )
    {
        // load Q into Mat3x3Packet
        AVX::Mat3x3Packet Q = AVX::load(
            b.q00, b.q01, b.q02,
            b.q10, b.q11, b.q12,
            b.q20, b.q21, b.q22
        );

        F.a00 = _mm256_sub_pd(
            _mm256_load_pd(b.p0x),
            _mm256_load_pd(b.p3x));

        F.a10 = _mm256_sub_pd(
            _mm256_load_pd(b.p0y),
            _mm256_load_pd(b.p3y));

        F.a20 = _mm256_sub_pd(
            _mm256_load_pd(b.p0z),
            _mm256_load_pd(b.p3z));
        
        F.a01 = _mm256_sub_pd(
            _mm256_load_pd(b.p1x),
            _mm256_load_pd(b.p3x));

        F.a11 = _mm256_sub_pd(
            _mm256_load_pd(b.p1y),
            _mm256_load_pd(b.p3y));

        F.a21 = _mm256_sub_pd(
            _mm256_load_pd(b.p1z),
            _mm256_load_pd(b.p3z));

        F.a02 = _mm256_sub_pd(
            _mm256_load_pd(b.p2x),
            _mm256_load_pd(b.p3x));

        F.a12 = _mm256_sub_pd(
            _mm256_load_pd(b.p2y),
            _mm256_load_pd(b.p3y));

        F.a22 = _mm256_sub_pd(
            _mm256_load_pd(b.p2z),
            _mm256_load_pd(b.p3z));

        AVX::matmul_right_inplace(F, Q);
    }

    static void accumulate4(
        unsigned e_idx[4],
        const NeoHookeanEnergyPool& energies,
        ParticlePool& particles,
        unsigned local_idx[4],
        Mat3r& particle_H,
        Vec3r& particle_G,
        Real dt
    )
    {
        // extract information into batch struct
        // which is SoA
        Batch4 batch;
        for (int lane = 0; lane < 4; lane++)
        {
            const auto& tet = energies.data[e_idx[lane]];
            unsigned li = local_idx[lane];

            Vec3r qi;

            if (li < 3)
            {
                qi = tet.Q.row(li);
                batch.sign[lane] = 1.0;
            }
            else
            {
                qi =
                    tet.Q.row(0)
                + tet.Q.row(1)
                + tet.Q.row(2);

                batch.sign[lane] = -1.0;
            }

            batch.qx[lane] = qi.x();
            batch.qy[lane] = qi.y();
            batch.qz[lane] = qi.z();

            batch.V[lane] = tet.rest_volume;
            batch.lambda[lane] = tet.lambda;
            batch.mu[lane] = tet.mu;
            batch.kd[lane] = tet.kd;

            const Vec3r& p0 = particles.positions[tet.particle_indices[0]];
            const Vec3r& p1 = particles.positions[tet.particle_indices[1]];
            const Vec3r& p2 = particles.positions[tet.particle_indices[2]];
            const Vec3r& p3 = particles.positions[tet.particle_indices[3]];
            batch.p0x[lane] = p0[0];
            batch.p0y[lane] = p0[1];
            batch.p0z[lane] = p0[2];

            batch.p1x[lane] = p1[0];
            batch.p1y[lane] = p1[1];
            batch.p1z[lane] = p1[2];

            batch.p2x[lane] = p2[0];
            batch.p2y[lane] = p2[1];
            batch.p2z[lane] = p2[2];

            batch.p3x[lane] = p3[0];
            batch.p3y[lane] = p3[1];
            batch.p3z[lane] = p3[2];

            batch.q00[lane] = tet.Q(0,0);
            batch.q01[lane] = tet.Q(0,1);
            batch.q02[lane] = tet.Q(0,2);
            batch.q10[lane] = tet.Q(1,0);
            batch.q11[lane] = tet.Q(1,1);
            batch.q12[lane] = tet.Q(1,2);
            batch.q20[lane] = tet.Q(2,0);
            batch.q21[lane] = tet.Q(2,1);
            batch.q22[lane] = tet.Q(2,2);

            batch.Eprev00[lane] = tet.E_prev(0,0);
            batch.Eprev01[lane] = tet.E_prev(0,1);
            batch.Eprev02[lane] = tet.E_prev(0,2);
            batch.Eprev10[lane] = tet.E_prev(1,0);
            batch.Eprev11[lane] = tet.E_prev(1,1);
            batch.Eprev12[lane] = tet.E_prev(1,2);
            batch.Eprev20[lane] = tet.E_prev(2,0);
            batch.Eprev21[lane] = tet.E_prev(2,1);
            batch.Eprev22[lane] = tet.E_prev(2,2);
        }

        AVX::Mat3x3Packet Eprev = AVX::load(
            batch.Eprev00, batch.Eprev01, batch.Eprev02,
            batch.Eprev10, batch.Eprev11, batch.Eprev12,
            batch.Eprev20, batch.Eprev21, batch.Eprev22
        );

        AVX::Mat3x3Packet F;
        computeF4(batch, F);

        AVX::Vec3Packet c01, c12, c20;
        AVX::cross_columns_all(F, c01, c12, c20);

        __m256d detF = AVX::determinant_packet(F, c12);


        AVX::Vec3Packet qi = AVX::load(batch.qx, batch.qy, batch.qz);
        AVX::Vec3Packet Fqi, Fcross_qi;
        AVX::matvec_packet(F, qi, Fqi);
        AVX::matvec_columns(c12, c20, c01, qi, Fcross_qi);
        
        __m256d V = _mm256_load_pd(batch.V);
        __m256d lambda = _mm256_load_pd(batch.lambda);
        __m256d mu = _mm256_load_pd(batch.mu);
        __m256d kd = _mm256_load_pd(batch.kd);
        __m256d sign = _mm256_load_pd(batch.sign);

        // hydrostatic multiplier
        __m256d hyd_mult = _mm256_mul_pd(sign, 
            _mm256_mul_pd(
                _mm256_mul_pd(V,lambda),
                _mm256_sub_pd(
                    _mm256_sub_pd(
                        detF,
                        _mm256_div_pd(mu,lambda)),
                    _mm256_set1_pd(1.0)))
        );

        // deviatoric multiplier
        __m256d dev_mult = _mm256_mul_pd( sign, _mm256_mul_pd(V, mu));

        // compute rate of change of Green strain
        AVX::Mat3x3Packet FFt;
        AVX::matmul_transpose_right(F, FFt);

        AVX::Mat3x3Packet FtF;
        AVX::matmul_transpose_left(F, FtF);
        
        __m256d qi_norm2 = AVX::squaredNorm3_packet(qi);

        // gradient
        __m256d damp_mult = _mm256_mul_pd( _mm256_set1_pd(2/dt), _mm256_mul_pd(V,kd) );
        __m256d signed_damp_mult = _mm256_mul_pd(sign, damp_mult);

        AVX::Vec3Packet G;

        // hyd*Fcross_qi + dev*Fqi
        G.x = _mm256_add_pd(
        _mm256_mul_pd(hyd_mult,Fcross_qi.x),
        _mm256_mul_pd(dev_mult,Fqi.x));

        G.y = _mm256_add_pd(
            _mm256_mul_pd(hyd_mult,Fcross_qi.y),
            _mm256_mul_pd(dev_mult,Fqi.y));

        G.z = _mm256_add_pd(
            _mm256_mul_pd(hyd_mult,Fcross_qi.z),
            _mm256_mul_pd(dev_mult,Fqi.z));

        // E_dot = ((1/2 *(FtF-I) -Eprev)/dt

        AVX::Mat3x3Packet Edot = FtF;
        const __m256d invdt =
            _mm256_set1_pd(1/dt);

        Edot.a00 = _mm256_mul_pd(
            _mm256_sub_pd(
                _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_sub_pd(Edot.a00,_mm256_set1_pd(1.0))),
                Eprev.a00),
            invdt);

        Edot.a01 = _mm256_mul_pd(
            _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), Edot.a01), Eprev.a01),
            invdt);

        Edot.a02 = _mm256_mul_pd(
            _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), Edot.a02), Eprev.a02),
            invdt);

        Edot.a10 = _mm256_mul_pd(
            _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), Edot.a10), Eprev.a10),
            invdt);

        Edot.a11 = _mm256_mul_pd(
            _mm256_sub_pd(
                _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_sub_pd(Edot.a11,_mm256_set1_pd(1.0))),
                Eprev.a11),
            invdt);

        Edot.a12 = _mm256_mul_pd(
            _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), Edot.a12), Eprev.a12),
            invdt);

        Edot.a20 = _mm256_mul_pd(
            _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), Edot.a20), Eprev.a20),
            invdt);

        Edot.a21 = _mm256_mul_pd(
            _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), Edot.a21), Eprev.a21),
            invdt);

        Edot.a22 = _mm256_mul_pd(
            _mm256_sub_pd(
                _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_sub_pd(Edot.a22,_mm256_set1_pd(1.0))),
                Eprev.a22),
            invdt);

        
        AVX::Vec3Packet Edot_qi,F_Edot_qi;
        AVX::matvec_packet(Edot, qi, Edot_qi);
        AVX::matvec_packet(F, Edot_qi, F_Edot_qi);

        G.x = _mm256_fmadd_pd(signed_damp_mult, F_Edot_qi.x, G.x);
        G.y = _mm256_fmadd_pd(signed_damp_mult, F_Edot_qi.y, G.y);
        G.z = _mm256_fmadd_pd(signed_damp_mult, F_Edot_qi.z, G.z);

        /** === Hessian === */


        AVX::Mat3x3Packet H;
        __m256d Vlam = _mm256_mul_pd( V, lambda );
        // hydrostatic
        AVX::scaled_outer_product(Fcross_qi, Vlam, H);

        // deviatoric + diagonal damping
        __m256d diag = _mm256_add_pd(
            _mm256_mul_pd(qi_norm2, _mm256_mul_pd(V, mu)),
            _mm256_mul_pd(damp_mult, AVX::dot3_packet(qi, Edot_qi))
        );
        AVX::add_diagonal(H, diag);

        // damping
        __m256d alpha = _mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(0.5), damp_mult), invdt);
        H.a00 = AVX::madd_outer(H.a00, alpha, qi_norm2, FFt.a00, Fqi.x, Fqi.x);
        H.a01 = AVX::madd_outer(H.a01, alpha, qi_norm2, FFt.a01, Fqi.x, Fqi.y);
        H.a02 = AVX::madd_outer(H.a02, alpha, qi_norm2, FFt.a02, Fqi.x, Fqi.z);

        H.a10 = AVX::madd_outer(H.a10, alpha, qi_norm2, FFt.a10, Fqi.y, Fqi.x);
        H.a11 = AVX::madd_outer(H.a11, alpha, qi_norm2, FFt.a11, Fqi.y, Fqi.y);
        H.a12 = AVX::madd_outer(H.a12, alpha, qi_norm2, FFt.a12, Fqi.y, Fqi.z);

        H.a20 = AVX::madd_outer(H.a20, alpha, qi_norm2, FFt.a20, Fqi.z, Fqi.x);
        H.a21 = AVX::madd_outer(H.a21, alpha, qi_norm2, FFt.a21, Fqi.z, Fqi.y);
        H.a22 = AVX::madd_outer(H.a22, alpha, qi_norm2, FFt.a22, Fqi.z, Fqi.z);

        // particle_H += (damp_mult_factor/dt) * (FFt * qi_norm2 + Fqi * Fqi.transpose());

        AVX::accumulate_Vector3Packet(particle_G, G);
        AVX::accumulate_Mat3x3Packet(particle_H, H);
    }

    /** Computes the Hessian and gradient for all of the particles affected by this constraint.
     * Updates the accumulated vertex Hessian and gradients.
     * @param c_idx : the energy index
     * @param energies : the memory pool for the energies
     * @param particles : the memory pool for the particles
     * @param vertex_Hs : a memory pool storing the per-vertex Hessians - this is updated by this function
     * @param vertex_Gs : a memory pool storing the per-vertex gradients - this is updated by this function
     * @param dt : the time step
     */
    // static void accumulateHessianGradient(
    //     unsigned c_idx,
    //     const NeoHookeanEnergyPool& energies,
    //     ParticlePool& particles,
    //     Mat3r* vertex_Hs,
    //     Vec3r* vertex_Gs,
    //     Real dt
    // )
    // {
    //     // extract data from energy pool
    //     const Vec4i& indices = energies.particle_indices[c_idx];
    //     Real V = energies.rest_volumes[c_idx];
    //     const Mat3r& Q = energies.Qs[c_idx];
    //     Real lambda = energies.lambdas[c_idx];
    //     Real mu = energies.mus[c_idx];
    //     Real kd = energies.kds[c_idx];

    //     // compute F for this timestep and the previous timestep
    //     Mat3r F = computeF(
    //         particles.positions[indices[0]],
    //         particles.positions[indices[1]],
    //         particles.positions[indices[2]],
    //         particles.positions[indices[3]],
    //         Q);
    //     Mat3r F_prev = computeF(
    //         particles.previous_positions[indices[0]],
    //         particles.previous_positions[indices[1]],
    //         particles.previous_positions[indices[2]],
    //         particles.previous_positions[indices[3]],
    //         Q);
        
    //     // compute rate of change of Green strain
    //     Mat3r E = F.transpose() * F - Mat3r::Identity();
    //     Mat3r E_prev = F_prev.transpose() * F_prev - Mat3r::Identity();
    //     Mat3r E_dot = 1/dt * (E - E_prev);

    //     // hydrostatic gradient
    //     Mat3r F_cross;
    //     F_cross.col(0) = F.col(1).cross(F.col(2));
    //     F_cross.col(1) = F.col(2).cross(F.col(0));
    //     F_cross.col(2) = F.col(0).cross(F.col(1));

    //     Mat3r hyd_grad_full = V*lambda * (F.determinant() - mu/lambda - 1) * F_cross * Q.transpose();

    //     // deviatoric gradient
    //     Mat3r dev_grad_full = V*mu * F * Q.transpose();

    //     // compute the Hessian and gradient for the first 3 particles
    //     for (int i = 0; i < 3; i++)
    //     {
    //         Vec3r qi = Q.row(i);
    //         Vec3r Fqi = F*qi;

    //         // gradient
    //         Vec3r hyd_grad_i = hyd_grad_full.col(i);
    //         Vec3r dev_grad_i = dev_grad_full.col(i);
    //         Vec3r damp_grad_i = 2/dt * V * kd * F * E_dot * qi;
    //         vertex_Gs[indices[i]] += hyd_grad_i + dev_grad_i + damp_grad_i;

    //         // Hessian
    //         Mat3r hyd_hess_i = V * lambda * detF_grad_full.col(i) * detF_grad_full.col(i).transpose();
    //         Mat3r dev_hess_i = V * mu * qi.squaredNorm() * Mat3r::Identity();
    //         Mat3r damp_hess_i = 2/dt * V * kd * ( (qi.transpose() * E_dot * qi) * Mat3r::Identity() + 1/dt * (Fqi * Fqi.transpose() + F*F.transpose() * qi.squaredNorm()) );
    //         vertex_Hs[indices[i]] += hyd_hess_i + dev_hess_i + damp_hess_i;
    //     }

    //     // compute Hessian and gradient for the 4th particle
    //     // gradient
    //     Vec3r q4 = Q.row(0) + Q.row(1) + Q.row(2);
    //     Vec3r hyd_grad_4 = -hyd_grad_full.col(0) - hyd_grad_full.col(1) - hyd_grad_full.col(2);
    //     Vec3r dev_grad_4 = -dev_grad_full.col(0) - dev_grad_full.col(1) - dev_grad_full.col(2);
    //     Vec3r damp_grad_4 = -2/dt * V * kd * F * E_dot * q4;
    //     vertex_Gs[indices[i]] += hyd_grad_4 + dev_grad_4 + damp_grad_4;

    //     // Hessian
    //     Vec3r a3 =
    //         -detF_grad_full.col(0)
    //         -detF_grad_full.col(1)
    //         -detF_grad_full.col(2);

    //     Mat3r hyd_hess =
    //         V*lambda * a3 * a3.transpose();

    //     Mat3r dev_hess =
    //         V*mu * q4.squaredNorm() * Mat3r::Identity();

    //     Vec3r Fq4 = F * q4;
    //     Mat3r damp_hess = 2/dt * V * kd * ( (q4.transpose() * E_dot * q4) * Mat3r::Identity() + 1/dt * (Fq4 * Fq4.transpose() + F*F.transpose() * q4.squaredNorm()) );
    //     vertex_Hs[indices[i]] += hyd_hess + dev_hess + damp_hess;
    // }
};

} // namespace Energy