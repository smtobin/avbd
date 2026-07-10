#pragma once

#include "common/common.hpp"
#include "immintrin.h"

struct AVX
{
struct alignas(32) Mat3x3Packet
{
    __m256d a00,a01,a02;
    __m256d a10,a11,a12;
    __m256d a20,a21,a22;
};


struct alignas(32) Vec3Packet
{
    __m256d x,y,z;
};

static inline Mat3x3Packet load(
    const double a00[4], const double a01[4], const double a02[4],
    const double a10[4], const double a11[4], const double a12[4],
    const double a20[4], const double a21[4], const double a22[4]
)
{
    Mat3x3Packet r;
    r.a00 = _mm256_load_pd(a00);
    r.a01 = _mm256_load_pd(a01);
    r.a02 = _mm256_load_pd(a02);
    r.a10 = _mm256_load_pd(a10);
    r.a11 = _mm256_load_pd(a11);
    r.a12 = _mm256_load_pd(a12);
    r.a20 = _mm256_load_pd(a20);
    r.a21 = _mm256_load_pd(a21);
    r.a22 = _mm256_load_pd(a22);
    return r;
}

static inline Vec3Packet load(
    const double x[4], const double y[4], const double z[4]
)
{
    Vec3Packet r;
    r.x = _mm256_load_pd(x);
    r.y = _mm256_load_pd(y);
    r.z = _mm256_load_pd(z);
    return r;
}

static inline double hsum_pd(__m256d x)
{
    __m128d hi = _mm256_extractf128_pd(x, 1);
    __m128d lo = _mm256_castpd256_pd128(x);

    lo = _mm_add_pd(lo, hi);
    lo = _mm_hadd_pd(lo, lo);

    return _mm_cvtsd_f64(lo);
}

static inline void accumulate_Vector3Packet(
    Vec3r& G,
    const Vec3Packet& Gp)
{
    G.x() += hsum_pd(Gp.x);
    G.y() += hsum_pd(Gp.y);
    G.z() += hsum_pd(Gp.z);
}

static inline void accumulate_Mat3x3Packet(
    Mat3r& H,
    const Mat3x3Packet& Hp)
{
    H(0,0) += hsum_pd(Hp.a00);
    H(0,1) += hsum_pd(Hp.a01);
    H(0,2) += hsum_pd(Hp.a02);
    H(1,0) += hsum_pd(Hp.a10);
    H(1,1) += hsum_pd(Hp.a11);
    H(1,2) += hsum_pd(Hp.a12);
    H(2,0) += hsum_pd(Hp.a20);
    H(2,1) += hsum_pd(Hp.a21);
    H(2,2) += hsum_pd(Hp.a22);
}

static inline __m256d fmadd3(
    __m256d a0, __m256d b0,
    __m256d a1, __m256d b1,
    __m256d a2, __m256d b2)
{
    // a0*b0 + a1*b1 + a2*b2
    return _mm256_fmadd_pd(
        a2, b2,
        _mm256_fmadd_pd(
            a1, b1,
            _mm256_mul_pd(a0,b0)));
}

static inline __m256d dot3_packet(
    const Vec3Packet& a,
    const Vec3Packet& b)
{
    return fmadd3(a.x, b.x, a.y, b.y, a.z, b.z);
}

static inline void matmul3x3_packet(
    const Mat3x3Packet& A,
    const Mat3x3Packet& B,
    Mat3x3Packet& C)
{
    C.a00 = fmadd3(
        A.a00,B.a00,
        A.a01,B.a10,
        A.a02,B.a20);

    C.a01 = fmadd3(
        A.a00,B.a01,
        A.a01,B.a11,
        A.a02,B.a21);

    C.a02 = fmadd3(
        A.a00,B.a02,
        A.a01,B.a12,
        A.a02,B.a22);



    C.a10 = fmadd3(
        A.a10,B.a00,
        A.a11,B.a10,
        A.a12,B.a20);

    C.a11 = fmadd3(
        A.a10,B.a01,
        A.a11,B.a11,
        A.a12,B.a21);

    C.a12 = fmadd3(
        A.a10,B.a02,
        A.a11,B.a12,
        A.a12,B.a22);



    C.a20 = fmadd3(
        A.a20,B.a00,
        A.a21,B.a10,
        A.a22,B.a20);

    C.a21 = fmadd3(
        A.a20,B.a01,
        A.a21,B.a11,
        A.a22,B.a21);

    C.a22 = fmadd3(
        A.a20,B.a02,
        A.a21,B.a12,
        A.a22,B.a22);
}

static inline void matmul_right_inplace(
    Mat3x3Packet& F,
    const Mat3x3Packet& Q)
{
    // Preserve original F values
    const __m256d f00 = F.a00;
    const __m256d f01 = F.a01;
    const __m256d f02 = F.a02;

    const __m256d f10 = F.a10;
    const __m256d f11 = F.a11;
    const __m256d f12 = F.a12;

    const __m256d f20 = F.a20;
    const __m256d f21 = F.a21;
    const __m256d f22 = F.a22;

    F.a00 = fmadd3(f00,Q.a00, f01,Q.a10, f02,Q.a20);
    F.a01 = fmadd3(f00,Q.a01, f01,Q.a11, f02,Q.a21);
    F.a02 = fmadd3(f00,Q.a02, f01,Q.a12, f02,Q.a22);

    F.a10 = fmadd3(f10,Q.a00, f11,Q.a10, f12,Q.a20);
    F.a11 = fmadd3(f10,Q.a01, f11,Q.a11, f12,Q.a21);
    F.a12 = fmadd3(f10,Q.a02, f11,Q.a12, f12,Q.a22);

    F.a20 = fmadd3(f20,Q.a00, f21,Q.a10, f22,Q.a20);
    F.a21 = fmadd3(f20,Q.a01, f21,Q.a11, f22,Q.a21);
    F.a22 = fmadd3(f20,Q.a02, f21,Q.a12, f22,Q.a22);
}

static inline void matmul_transpose_right(
    const Mat3x3Packet& F,
    Mat3x3Packet& FFt)
{
    // row 0 dot row 0
    FFt.a00 = fmadd3(
        F.a00,F.a01,F.a02,
        F.a00,F.a01,F.a02);

    // row 0 dot row 1
    FFt.a01 = fmadd3(
        F.a00,F.a01,F.a02,
        F.a10,F.a11,F.a12);

    // row 0 dot row 2
    FFt.a02 = fmadd3(
        F.a00,F.a01,F.a02,
        F.a20,F.a21,F.a22);


    // row 1 dot row 1
    FFt.a11 = fmadd3(
        F.a10,F.a11,F.a12,
        F.a10,F.a11,F.a12);

    // row 1 dot row 2
    FFt.a12 = fmadd3(
        F.a10,F.a11,F.a12,
        F.a20,F.a21,F.a22);


    // row 2 dot row 2
    FFt.a22 = fmadd3(
        F.a20,F.a21,F.a22,
        F.a20,F.a21,F.a22);


    // Symmetry
    FFt.a10 = FFt.a01;
    FFt.a20 = FFt.a02;
    FFt.a21 = FFt.a12;
}

static inline void matvec_columns(
    const Vec3Packet& c0,
    const Vec3Packet& c1,
    const Vec3Packet& c2,
    const Vec3Packet& q,
    Vec3Packet& r)
{
    r.x = _mm256_fmadd_pd(
        c2.x, q.z,
        _mm256_fmadd_pd(
            c1.x, q.y,
            _mm256_mul_pd(c0.x, q.x)));

    r.y = _mm256_fmadd_pd(
        c2.y, q.z,
        _mm256_fmadd_pd(
            c1.y, q.y,
            _mm256_mul_pd(c0.y, q.x)));

    r.z = _mm256_fmadd_pd(
        c2.z, q.z,
        _mm256_fmadd_pd(
            c1.z, q.y,
            _mm256_mul_pd(c0.z, q.x)));
}

static inline void matvec_packet(
    const Mat3x3Packet& F,
    const Vec3Packet& q,
    Vec3Packet& r)
{

    r.x = _mm256_fmadd_pd(
        F.a02, q.z,
        _mm256_fmadd_pd(
            F.a01, q.y,
            _mm256_mul_pd(F.a00,q.x)));

    r.y = _mm256_fmadd_pd(
        F.a12, q.z,
        _mm256_fmadd_pd(
            F.a11, q.y,
            _mm256_mul_pd(F.a10,q.x)));

    r.z = _mm256_fmadd_pd(
        F.a22, q.z,
        _mm256_fmadd_pd(
            F.a21, q.y,
            _mm256_mul_pd(F.a20,q.x)));
}

static inline void cross3_packet(
    const Vec3Packet& a,
    const Vec3Packet& b,
    Vec3Packet& c)
{
    c.x = _mm256_fnmadd_pd(
        a.z,
        b.y,
        _mm256_mul_pd(a.y,b.z));

    c.y = _mm256_fnmadd_pd(
        a.x,
        b.z,
        _mm256_mul_pd(a.z,b.x));

    c.z = _mm256_fnmadd_pd(
        a.y,
        b.x,
        _mm256_mul_pd(a.x,b.y));
}

static inline void cross_columns_1_2(
    const Mat3x3Packet& F,
    Vec3Packet& c)
{
    c.x = _mm256_fnmadd_pd(
        F.a21,
        F.a12,
        _mm256_mul_pd(F.a11,F.a22));

    c.y = _mm256_fnmadd_pd(
        F.a01,
        F.a22,
        _mm256_mul_pd(F.a21,F.a02));

    c.z = _mm256_fnmadd_pd(
        F.a11,
        F.a02,
        _mm256_mul_pd(F.a01,F.a12));
}

static inline void cross_columns_2_0(
    const Mat3x3Packet& F,
    Vec3Packet& c)
{
    c.x = _mm256_fnmadd_pd(
        F.a22,
        F.a10,
        _mm256_mul_pd(F.a12,F.a20));

    c.y = _mm256_fnmadd_pd(
        F.a02,
        F.a20,
        _mm256_mul_pd(F.a22,F.a00));

    c.z = _mm256_fnmadd_pd(
        F.a12,
        F.a00,
        _mm256_mul_pd(F.a02,F.a10));
}

static inline void cross_columns_0_1(
    const Mat3x3Packet& F,
    Vec3Packet& c)
{
    c.x = _mm256_fnmadd_pd(
        F.a20,
        F.a11,
        _mm256_mul_pd(F.a10,F.a21));

    c.y = _mm256_fnmadd_pd(
        F.a00,
        F.a21,
        _mm256_mul_pd(F.a20,F.a01));

    c.z = _mm256_fnmadd_pd(
        F.a10,
        F.a01,
        _mm256_mul_pd(F.a00,F.a11));
}

static inline void cross_columns_all(
    const Mat3x3Packet& F,
    Vec3Packet& c01,
    Vec3Packet& c12,
    Vec3Packet& c20)
{
    // Load columns
    const __m256d a00 = F.a00;
    const __m256d a01 = F.a01;
    const __m256d a02 = F.a02;

    const __m256d a10 = F.a10;
    const __m256d a11 = F.a11;
    const __m256d a12 = F.a12;

    const __m256d a20 = F.a20;
    const __m256d a21 = F.a21;
    const __m256d a22 = F.a22;


    // col0 x col1
    c01.x = _mm256_fnmadd_pd(
        a20, a11,
        _mm256_mul_pd(a10, a21));

    c01.y = _mm256_fnmadd_pd(
        a00, a21,
        _mm256_mul_pd(a20, a01));

    c01.z = _mm256_fnmadd_pd(
        a10, a01,
        _mm256_mul_pd(a00, a11));


    // col1 x col2
    c12.x = _mm256_fnmadd_pd(
        a21, a12,
        _mm256_mul_pd(a11, a22));

    c12.y = _mm256_fnmadd_pd(
        a01, a22,
        _mm256_mul_pd(a21, a02));

    c12.z = _mm256_fnmadd_pd(
        a11, a02,
        _mm256_mul_pd(a01, a12));


    // col2 x col0
    c20.x = _mm256_fnmadd_pd(
        a22, a10,
        _mm256_mul_pd(a12, a20));

    c20.y = _mm256_fnmadd_pd(
        a02, a20,
        _mm256_mul_pd(a22, a00));

    c20.z = _mm256_fnmadd_pd(
        a12, a00,
        _mm256_mul_pd(a02, a10));
}

static inline __m256d determinant_packet(
    const Mat3x3Packet& F,
    const Vec3Packet& c01)
{
    return _mm256_fmadd_pd(
        F.a20, c01.z,
        _mm256_fmadd_pd(
            F.a10, c01.y,
            _mm256_mul_pd(F.a00, c01.x)));
}

static inline __m256d squaredNorm3_packet(
    const Vec3Packet& q)
{
    return _mm256_fmadd_pd(
        q.z, q.z,
        _mm256_fmadd_pd(
            q.y, q.y,
            _mm256_mul_pd(q.x,q.x)));
}

static inline void scaled_outer_product(
    const Vec3Packet& v,
    __m256d scale,
    Mat3x3Packet& H)
{
    H.a00 = _mm256_mul_pd(
        scale,
        _mm256_mul_pd(v.x,v.x));

    H.a01 = _mm256_mul_pd(
        scale,
        _mm256_mul_pd(v.x,v.y));

    H.a02 = _mm256_mul_pd(
        scale,
        _mm256_mul_pd(v.x,v.z));

    H.a11 = _mm256_mul_pd(
        scale,
        _mm256_mul_pd(v.y,v.y));

    H.a12 = _mm256_mul_pd(
        scale,
        _mm256_mul_pd(v.y,v.z));

    H.a22 = _mm256_mul_pd(
        scale,
        _mm256_mul_pd(v.z,v.z));
    
    // exploit symmetry
    H.a10 = H.a01;
    H.a20 = H.a02;
    H.a21 = H.a12;
}

static inline void add_diagonal(
    Mat3x3Packet& H,
    __m256d d)
{
    H.a00 = _mm256_add_pd(H.a00, d);
    H.a11 = _mm256_add_pd(H.a11, d);
    H.a22 = _mm256_add_pd(H.a22, d);
}

static inline __m256d madd_outer(
    __m256d Hij,
    __m256d alpha,
    __m256d q2,
    __m256d Aij,
    __m256d ui,
    __m256d uj)
{
    // Hij += alpha * (q2*Aij + ui*uj)

    __m256d t =
        _mm256_fmadd_pd(
            q2,
            Aij,
            _mm256_mul_pd(ui,uj));

    return _mm256_fmadd_pd(alpha, t, Hij);
}

};

