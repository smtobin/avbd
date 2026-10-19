
#pragma once

#include <cmath>
#include <cassert>

#ifdef __CUDACC__
    #define HD __host__ __device__
#else
    #define HD
#endif

/** Simple Vec3 struct for GPU operations.
 * TODO: (10/19/26) Should the CPU code us this math layer as well?? (Probably yes long term)
 */
struct DVec3
{
    float x, y, z;

    // constructors
    HD constexpr DVec3() : x(0.0f), y(0.0f), z(0.0f) {}

    HD constexpr explicit DVec3(float s) : x(s), y(s), z(s) {}

    HD constexpr DVec3(float x_, float y_, float z_)
        : x(x_), y(y_), z(z_) {}

    // indexed access
    HD constexpr float& operator[](int i)
    {
        assert(i >= 0 && i < 3);
        return (&x)[i];
    }

    HD constexpr const float& operator[](int i) const
    {
        assert(i >= 0 && i < 3);
        return (&x)[i];
    }

    // unary operators
    HD constexpr DVec3 operator+() const { return *this; }
    HD constexpr DVec3 operator-() const { return {-x, -y, -z}; }

    // vector arithmetic
    HD constexpr DVec3 operator+(const DVec3& b) const
    {
        return {x + b.x, y + b.y, z + b.z};
    }

    HD constexpr DVec3 operator-(const DVec3& b) const
    {
        return {x - b.x, y - b.y, z - b.z};
    }

    HD constexpr DVec3 operator*(float s) const
    {
        return {x * s, y * s, z * s};
    }

    HD constexpr DVec3 operator/(float s) const
    {
        return {x / s, y / s, z / s};
    }

    // compound assignment
    HD constexpr DVec3& operator+=(const DVec3& b)
    {
        x += b.x; y += b.y; z += b.z;
        return *this;
    }

    HD constexpr DVec3& operator-=(const DVec3& b)
    {
        x -= b.x; y -= b.y; z -= b.z;
        return *this;
    }

    HD constexpr DVec3& operator*=(float s)
    {
        x *= s; y *= s; z *= s;
        return *this;
    }

    HD constexpr DVec3& operator/=(float s)
    {
        x /= s; y /= s; z /= s;
        return *this;
    }

    // geometric operations
    HD constexpr float dot(const DVec3& b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }

    HD constexpr DVec3 cross(const DVec3& b) const
    {
        return {
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x
        };
    }

    HD constexpr float squaredNorm() const
    {
        return dot(*this);
    }

    HD float norm() const
    {
        return ::sqrtf(squaredNorm());
    }

    // normalize this vector in place.
    // returns false if the vector is zero or non-finite.
    HD bool normalize()
    {
        const float n2 = squaredNorm();

        if (!(n2 > 0.0f) || !::isfinite(n2))
            return false;

        const float invN = 1.0f / ::sqrtf(n2);
        x *= invN;
        y *= invN;
        z *= invN;
        return true;
    }

    // return a normalized copy.
    // returns zero if normalization fails.
    HD DVec3 normalized() const
    {
        DVec3 result = *this;
        if (!result.normalize())
            return DVec3{};
        return result;
    }

    // component-wise multiplication and division
    HD constexpr DVec3 cwiseProduct(const DVec3& b) const
    {
        return {x * b.x, y * b.y, z * b.z};
    }

    HD constexpr DVec3 cwiseQuotient(const DVec3& b) const
    {
        return {x / b.x, y / b.y, z / b.z};
    }

    // Component-wise min/max
    HD DVec3 cwiseMin(const DVec3& b) const
    {
        return {
            ::fminf(x, b.x),
            ::fminf(y, b.y),
            ::fminf(z, b.z)
        };
    }

    HD DVec3 cwiseMax(const DVec3& b) const
    {
        return {
            ::fmaxf(x, b.x),
            ::fmaxf(y, b.y),
            ::fmaxf(z, b.z)
        };
    }

    HD constexpr float maxCoeff() const
    {
        return x > y ? (x > z ? x : z) : (y > z ? y : z);
    }

    HD constexpr float minCoeff() const
    {
        return x < y ? (x < z ? x : z) : (y < z ? y : z);
    }
};

// scalar multiplication in either order
HD constexpr inline DVec3 operator*(float s, const DVec3& v)
{
    return v * s;
}
