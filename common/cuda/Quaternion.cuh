
#pragma once

#include "Vec3.cuh"
#include <cmath>
#include <cassert>

#ifdef __CUDACC__
    #define HD __host__ __device__
#else
    #define HD
#endif

/** Quaternion class for GPU operations.
 * TODO: (10/19/26) CPU code should use this struct as well
 */
struct DQuaternion
{
    // Storage order: x, y, z, w
    float x, y, z, w;

    // Constructors
    HD constexpr DQuaternion()
        : x(0.0f), y(0.0f), z(0.0f), w(1.0f) {}

    // Eigen-style constructor: (w, x, y, z)
    HD constexpr DQuaternion(float w_, float x_, float y_, float z_)
        : x(x_), y(y_), z(z_), w(w_) {}

    HD constexpr DQuaternion(const DVec3& vec, float scalar)
        : x(vec.x), y(vec.y), z(vec.z), w(scalar) {}

    // Identity DQuaternion
    HD static constexpr DQuaternion Identity()
    {
        return DQuaternion(1.0f, 0.0f, 0.0f, 0.0f);
    }

    // Access vector/scalar parts
    HD constexpr DVec3 vec() const { return {x, y, z}; }
    HD constexpr float scalar() const { return w; }

    // Unary operations
    HD constexpr DQuaternion operator+() const { return *this; }

    HD constexpr DQuaternion operator-() const
    {
        return DQuaternion(-w, -x, -y, -z);
    }

    // Addition and subtraction are coefficient-wise
    HD constexpr DQuaternion operator+(const DQuaternion& b) const
    {
        return DQuaternion(w + b.w, x + b.x, y + b.y, z + b.z);
    }

    HD constexpr DQuaternion operator-(const DQuaternion& b) const
    {
        return DQuaternion(w - b.w, x - b.x, y - b.y, z - b.z);
    }

    // Scalar multiplication/division
    HD constexpr DQuaternion operator*(float s) const
    {
        return DQuaternion(w * s, x * s, y * s, z * s);
    }

    HD constexpr DQuaternion operator/(float s) const
    {
        return DQuaternion(w / s, x / s, y / s, z / s);
    }

    HD constexpr DQuaternion& operator+=(const DQuaternion& b)
    {
        x += b.x; y += b.y; z += b.z; w += b.w;
        return *this;
    }

    HD constexpr DQuaternion& operator-=(const DQuaternion& b)
    {
        x -= b.x; y -= b.y; z -= b.z; w -= b.w;
        return *this;
    }

    HD constexpr DQuaternion& operator*=(float s)
    {
        x *= s; y *= s; z *= s; w *= s;
        return *this;
    }

    HD constexpr DQuaternion& operator/=(float s)
    {
        x /= s; y /= s; z /= s; w /= s;
        return *this;
    }

    // DQuaternion product (Hamilton product)
    HD constexpr DQuaternion operator*(const DQuaternion& b) const
    {
        return DQuaternion(
            w*b.w - x*b.x - y*b.y - z*b.z,
            w*b.x + x*b.w + y*b.z - z*b.y,
            w*b.y - x*b.z + y*b.w + z*b.x,
            w*b.z + x*b.y - y*b.x + z*b.w
        );
    }

    // Dot product and norms
    HD constexpr float dot(const DQuaternion& b) const
    {
        return x*b.x + y*b.y + z*b.z + w*b.w;
    }

    HD constexpr float squaredNorm() const
    {
        return dot(*this);
    }

    HD float norm() const
    {
        return ::sqrtf(squaredNorm());
    }

    // Normalize in place; returns false for zero or non-finite norm.
    HD bool normalize()
    {
        const float n2 = squaredNorm();

        if (!(n2 > 0.0f) || !::isfinite(n2))
            return false;

        const float invN = 1.0f / ::sqrtf(n2);
        x *= invN;
        y *= invN;
        z *= invN;
        w *= invN;
        return true;
    }

    HD DQuaternion normalized() const
    {
        DQuaternion q = *this;
        if (!q.normalize())
            return Identity();
        return q;
    }

    // Conjugate: inverse for unit DQuaternions
    HD constexpr DQuaternion conjugate() const
    {
        return DQuaternion(w, -x, -y, -z);
    }

    // General inverse: conjugate(q) / squaredNorm(q)
    HD DQuaternion inverse() const
    {
        const float n2 = squaredNorm();

        if (!(n2 > 0.0f) || !::isfinite(n2))
            return Identity();

        return conjugate() / n2;
    }

    // Rotate a vector. This formula assumes q is unit length.
    HD DVec3 operator*(const DVec3& v) const
    {
        const DVec3 qv{x, y, z};
        const DVec3 t = 2.0f * qv.cross(v);
        return v + w * t + qv.cross(t);
    }

    // Angle in radians, assuming q is normalized.
    HD float angularDistance(const DQuaternion& b) const
    {
        // q and -q represent the same rotation.
        float d = ::fabsf(dot(b));
        d = d < 1.0f ? d : 1.0f;
        return 2.0f * ::acosf(d);
    }

    // Construct a unit DQuaternion from an axis and angle in radians.
    // Returns identity if the axis cannot be normalized.
    HD static DQuaternion FromAxisAngle(DVec3 axis, float angle)
    {
        if (!axis.normalize())
            return Identity();

        const float half = 0.5f * angle;
        const float s = ::sinf(half);
        const float c = ::cosf(half);

        return DQuaternion(c, axis.x * s, axis.y * s, axis.z * s);
    }

    // Convert to a 3x3 rotation matrix, row-major output.
    // q must be unit length.
    HD void toRotationMatrix(float out[9]) const
    {
        const float xx = x*x, yy = y*y, zz = z*z;
        const float xy = x*y, xz = x*z, yz = y*z;
        const float wx = w*x, wy = w*y, wz = w*z;

        out[0] = 1.0f - 2.0f * (yy + zz);
        out[1] = 2.0f * (xy - wz);
        out[2] = 2.0f * (xz + wy);

        out[3] = 2.0f * (xy + wz);
        out[4] = 1.0f - 2.0f * (xx + zz);
        out[5] = 2.0f * (yz - wx);

        out[6] = 2.0f * (xz - wy);
        out[7] = 2.0f * (yz + wx);
        out[8] = 1.0f - 2.0f * (xx + yy);
    }
};

// Scalar multiplication in either order
HD constexpr inline DQuaternion operator*(float s, const DQuaternion& q)
{
    return q * s;
}
