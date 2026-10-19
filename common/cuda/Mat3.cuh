
#pragma once

#include "common/cuda/Vec3.cuh"

#include <cmath>
#include <cstddef>

#ifndef HD
    #ifdef __CUDACC__
        #define HD __host__ __device__
    #else
        #define HD
    #endif
#endif

struct DMat3
{
    using Scalar = float;

    // Row-major storage:
    // [ m[0] m[1] m[2] ]
    // [ m[3] m[4] m[5] ]
    // [ m[6] m[7] m[8] ]
    Scalar m[9];

    // Constructors
    HD DMat3() = default;

    HD explicit DMat3(Scalar diagonal)
        : m{
            diagonal, 0.0f, 0.0f,
            0.0f, diagonal, 0.0f,
            0.0f, 0.0f, diagonal
        }
    {}

    HD DMat3(
        Scalar m00, Scalar m01, Scalar m02,
        Scalar m10, Scalar m11, Scalar m12,
        Scalar m20, Scalar m21, Scalar m22)
        : m{
            m00, m01, m02,
            m10, m11, m12,
            m20, m21, m22
        }
    {}

    // Factories
    HD static DMat3 Identity() {
        return DMat3(1.0f);
    }

    HD static DMat3 Zero() {
        return DMat3(0.0f);
    }

    HD static DMat3 Diagonal(const DVec3& d) {
        return DMat3(
            d.x, 0.0f, 0.0f,
            0.0f, d.y, 0.0f,
            0.0f, 0.0f, d.z
        );
    }

    // Construct a matrix from row or column vectors.
    HD static DMat3 FromRows(
        const DVec3& r0,
        const DVec3& r1,
        const DVec3& r2)
    {
        return DMat3(
            r0.x, r0.y, r0.z,
            r1.x, r1.y, r1.z,
            r2.x, r2.y, r2.z
        );
    }

    HD static DMat3 FromColumns(
        const DVec3& c0,
        const DVec3& c1,
        const DVec3& c2)
    {
        return DMat3(
            c0.x, c1.x, c2.x,
            c0.y, c1.y, c2.y,
            c0.z, c1.z, c2.z
        );
    }

    // Element access (zero-based, unchecked).
    HD Scalar& operator()(int row, int col) {
        return m[row * 3 + col];
    }

    HD const Scalar& operator()(int row, int col) const {
        return m[row * 3 + col];
    }

    HD Scalar& operator[](int i) {
        return m[i];
    }

    HD const Scalar& operator[](int i) const {
        return m[i];
    }

    // Row and column access.
    HD DVec3 row(int i) const {
        return DVec3(m[i * 3], m[i * 3 + 1], m[i * 3 + 2]);
    }

    HD DVec3 col(int i) const {
        return DVec3(m[i], m[i + 3], m[i + 6]);
    }

    HD void setRow(int i, const DVec3& v) {
        m[i * 3]     = v.x;
        m[i * 3 + 1] = v.y;
        m[i * 3 + 2] = v.z;
    }

    HD void setCol(int i, const DVec3& v) {
        m[i]     = v.x;
        m[i + 3] = v.y;
        m[i + 6] = v.z;
    }

    // Unary operations.
    HD DMat3 operator+() const {
        return *this;
    }

    HD DMat3 operator-() const {
        DMat3 result;
        for (int i = 0; i < 9; ++i)
            result.m[i] = -m[i];
        return result;
    }

    // Matrix arithmetic.
    HD DMat3 operator+(const DMat3& b) const {
        DMat3 result;
        for (int i = 0; i < 9; ++i)
            result.m[i] = m[i] + b.m[i];
        return result;
    }

    HD DMat3 operator-(const DMat3& b) const {
        DMat3 result;
        for (int i = 0; i < 9; ++i)
            result.m[i] = m[i] - b.m[i];
        return result;
    }

    HD DMat3& operator+=(const DMat3& b) {
        for (int i = 0; i < 9; ++i)
            m[i] += b.m[i];
        return *this;
    }

    HD DMat3& operator-=(const DMat3& b) {
        for (int i = 0; i < 9; ++i)
            m[i] -= b.m[i];
        return *this;
    }

    // Scalar arithmetic.
    HD DMat3 operator*(Scalar s) const {
        DMat3 result;
        for (int i = 0; i < 9; ++i)
            result.m[i] = m[i] * s;
        return result;
    }

    HD DMat3 operator/(Scalar s) const {
        const Scalar inv = 1.0f / s;
        return *this * inv;
    }

    HD DMat3& operator*=(Scalar s) {
        for (int i = 0; i < 9; ++i)
            m[i] *= s;
        return *this;
    }

    HD DMat3& operator/=(Scalar s) {
        const Scalar inv = 1.0f / s;
        return *this *= inv;
    }

    // Matrix-vector multiplication: column-vector convention.
    HD DVec3 operator*(const DVec3& v) const {
        return DVec3(
            m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[3] * v.x + m[4] * v.y + m[5] * v.z,
            m[6] * v.x + m[7] * v.y + m[8] * v.z
        );
    }

    // Matrix-matrix multiplication.
    HD DMat3 operator*(const DMat3& b) const {
        DMat3 result = Zero();

        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                for (int k = 0; k < 3; ++k) {
                    result(r, c) += (*this)(r, k) * b(k, c);
                }
            }
        }

        return result;
    }

    HD DMat3& operator*=(const DMat3& b) {
        *this = *this * b;
        return *this;
    }

    // Transpose.
    HD DMat3 transpose() const {
        return DMat3(
            m[0], m[3], m[6],
            m[1], m[4], m[7],
            m[2], m[5], m[8]
        );
    }

    // Basic matrix properties.
    HD Scalar trace() const {
        return m[0] + m[4] + m[8];
    }

    HD Scalar determinant() const {
        return
            m[0] * (m[4] * m[8] - m[5] * m[7])
          - m[1] * (m[3] * m[8] - m[5] * m[6])
          + m[2] * (m[3] * m[7] - m[4] * m[6]);
    }

    HD Scalar squaredNorm() const {
        Scalar result = 0.0f;
        for (int i = 0; i < 9; ++i)
            result += m[i] * m[i];
        return result;
    }

    HD Scalar norm() const {
        return sqrtf(squaredNorm());
    }

    // Inverse using the adjugate.
    //
    // Returns false if the determinant is too close to zero.
    // On failure, *this is unchanged.
    HD bool tryInverse(DMat3& result, Scalar eps = 1e-8f) const {
        const Scalar det = determinant();

        if (!isfinite(det) || fabsf(det) <= eps)
            return false;

        const Scalar invDet = 1.0f / det;

        DMat3 inverse(
             m[4] * m[8] - m[5] * m[7],
             m[2] * m[7] - m[1] * m[8],
             m[1] * m[5] - m[2] * m[4],

             m[5] * m[6] - m[3] * m[8],
             m[0] * m[8] - m[2] * m[6],
             m[2] * m[3] - m[0] * m[5],

             m[3] * m[7] - m[4] * m[6],
             m[1] * m[6] - m[0] * m[7],
             m[0] * m[4] - m[1] * m[3]
        );

        result = inverse * invDet;
        return true;
    }

    // Skew-symmetric cross-product matrix:
    // crossMatrix(v) * w == v.cross(w)
    HD static DMat3 CrossMatrix(const DVec3& v) {
        return DMat3(
             0.0f, -v.z,  v.y,
              v.z,  0.0f, -v.x,
             -v.y,  v.x,  0.0f
        );
    }

    // Outer product: a * b^T.
    HD static DMat3 OuterProduct(const DVec3& a, const DVec3& b) {
        return DMat3(
            a.x * b.x, a.x * b.y, a.x * b.z,
            a.y * b.x, a.y * b.y, a.y * b.z,
            a.z * b.x, a.z * b.y, a.z * b.z
        );
    }
};

// Scalar on the left: s * M.
HD inline DMat3 operator*(DMat3::Scalar s, const DMat3& a) {
    return a * s;
}
