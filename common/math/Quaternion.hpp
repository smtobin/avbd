
#pragma once

#include "common/math/MatN.hpp"
#include <cmath>
#include <cassert>

#ifdef __CUDACC__
    #define HD __host__ __device__
#else
    #define HD
#endif

/** QuaternionT class for GPU operations.
 * TODO: (10/19/26) CPU code should use this struct as _well
 */
template <typename Real = float>
class QuaternionT
{
private:
    using Vec3r = VecN<3, Real>;
    // Storage order: _x, _y, _z, _w

    Real _x, _y, _z, _w;
public:
    // Constructors
    HD constexpr QuaternionT()
        : _x(0.0f), _y(0.0f), _z(0.0f), _w(1.0f) {}

    // Eigen-style constructor: (w, _x, _y, _z)
    HD constexpr QuaternionT(Real w, Real x, Real y, Real z)
        : _x(x), _y(y), _z(z), _w(w) {}

    HD constexpr QuaternionT(const Vec3r& vec, Real scalar)
        : _x(vec.x), _y(vec.y), _z(vec.z), _w(scalar) {}

    // Identity QuaternionT
    HD static constexpr QuaternionT Identity()
    {
        return QuaternionT(1.0f, 0.0f, 0.0f, 0.0f);
    }

    // Access vector/scalar parts
    HD constexpr Vec3r vec() const { return {_x, _y, _z}; }
    HD constexpr Real scalar() const { return _w; }

    HD Real& x() { return _x; }
    HD const Real& x() const { return _x; }

    HD Real& y() { return _y; }
    HD const Real& y() const { return _y; }

    HD Real& z() { return _z; }
    HD const Real& z() const { return _z; }

    HD Real& w() { return _w; }
    HD const Real& w() const { return _w; }

    // Unary operations
    HD constexpr QuaternionT operator+() const { return *this; }

    HD constexpr QuaternionT operator-() const
    {
        return QuaternionT(-w, -x, -y, -z);
    }

    // Addition and subtraction are coefficient-wise
    HD constexpr QuaternionT operator+(const QuaternionT& b) const
    {
        return QuaternionT(_w + b.w, _x + b.x, _y + b.y, _z + b.z);
    }

    HD constexpr QuaternionT operator-(const QuaternionT& b) const
    {
        return QuaternionT(_w - b.w, _x - b.x, _y - b.y, _z - b.z);
    }

    // Scalar multiplication/division
    HD constexpr QuaternionT operator*(Real s) const
    {
        return QuaternionT(_w * s, _x * s, _y * s, _z * s);
    }

    HD constexpr QuaternionT operator/(Real s) const
    {
        return QuaternionT(_w / s, _x / s, _y / s, _z / s);
    }

    HD constexpr QuaternionT& operator+=(const QuaternionT& b)
    {
        _x += b.x; _y += b.y; _z += b.z; _w += b.w;
        return *this;
    }

    HD constexpr QuaternionT& operator-=(const QuaternionT& b)
    {
        _x -= b.x; _y -= b.y; _z -= b.z; _w -= b.w;
        return *this;
    }

    HD constexpr QuaternionT& operator*=(Real s)
    {
        _x *= s; _y *= s; _z *= s; _w *= s;
        return *this;
    }

    HD constexpr QuaternionT& operator/=(Real s)
    {
        _x /= s; _y /= s; _z /= s; _w /= s;
        return *this;
    }

    // QuaternionT product (Hamilton product)
    HD constexpr QuaternionT operator*(const QuaternionT& b) const
    {
        return QuaternionT(
            _w*b._w - _x*b._x - _y*b._y - _z*b._z,
            _w*b._x + _x*b._w + _y*b._z - _z*b._y,
            _w*b._y - _x*b._z + _y*b._w + _z*b._x,
            _w*b._z + _x*b._y - _y*b._x + _z*b._w
        );
    }

    // Dot product and norms
    HD constexpr Real dot(const QuaternionT& b) const
    {
        return _x*b._x + _y*b._y + _z*b._z + _w*b._w;
    }

    HD constexpr Real squaredNorm() const
    {
        return dot(*this);
    }

    HD Real norm() const
    {
        return ::sqrtf(squaredNorm());
    }

    // Normalize in place; returns false for _zero or non-finite norm.
    HD bool normalize()
    {
        const Real n2 = squaredNorm();

        if (!(n2 > 0.0f) || !::isfinite(n2))
            return false;

        const Real invN = 1.0f / ::sqrtf(n2);
        _x *= invN;
        _y *= invN;
        _z *= invN;
        _w *= invN;
        return true;
    }

    HD QuaternionT normalized() const
    {
        QuaternionT q = *this;
        if (!q.normalize())
            return Identity();
        return q;
    }

    // Conjugate: inverse for unit QuaternionTs
    HD constexpr QuaternionT conjugate() const
    {
        return QuaternionT(_w, -_x, -_y, -_z);
    }

    // General inverse: conjugate(q) / squaredNorm(q)
    HD QuaternionT inverse() const
    {
        const Real n2 = squaredNorm();

        if (!(n2 > 0.0f) || !::isfinite(n2))
            return Identity();

        return conjugate() / n2;
    }

    // Rotate a vector. This formula assumes q is unit length.
    HD Vec3r operator*(const Vec3r& v) const
    {
        const Vec3r qv{_x, _y, _z};
        const Vec3r t = 2.0f * qv.cross(v);
        return v + _w * t + qv.cross(t);
    }

    // Angle in radians, assuming q is normalized.
    HD Real angularDistance(const QuaternionT& b) const
    {
        // q and -q represent the same rotation.
        Real d = ::fabsf(dot(b));
        d = d < 1.0f ? d : 1.0f;
        return 2.0f * ::acosf(d);
    }

    // Construct a unit QuaternionT from an axis and angle in radians.
    // Returns identity if the axis cannot be normalized.
    HD static QuaternionT FromAxisAngle(Vec3r axis, Real angle)
    {
        if (!axis.normalize())
            return Identity();

        const Real half = 0.5f * angle;
        const Real s = ::sinf(half);
        const Real c = ::cosf(half);

        return QuaternionT(c, axis.x * s, axis.y * s, axis.z * s);
    }

    // Convert to a 3x3 rotation matrix, row-major output.
    // q must be unit length.
    HD Mat3<Real> toRotationMatrix() const
    {
        const Real _xx = _x*_x, _yy = _y*_y, _zz = _z*_z;
        const Real _xy = _x*_y, _xz = _x*_z, _yz = _y*_z;
        const Real _wx = _w*_x, _wy = _w*_y, _wz = _w*_z;

        Mat3<Real> out;
        out(0,0) = 1.0f - 2.0f * (_yy + _zz);
        out(0,1) = 2.0f * (_xy - _wz);
        out(0,2) = 2.0f * (_xz + _wy);

        out(1,0) = 2.0f * (_xy + _wz);
        out(1,1) = 1.0f - 2.0f * (_xx + _zz);
        out(1,2) = 2.0f * (_yz - _wx);

        out(2,0) = 2.0f * (_xz - _wy);
        out(2,1) = 2.0f * (_yz + _wx);
        out(2,2) = 1.0f - 2.0f * (_xx + _yy);

        return out;
    }
};

// Scalar multiplication in either order
template <typename Real>
HD constexpr inline QuaternionT<Real> operator*(Real s, const QuaternionT<Real>& q)
{
    return q * s;
}
