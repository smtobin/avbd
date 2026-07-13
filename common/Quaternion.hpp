#pragma once

#include "common/common.hpp"

struct Quaternion
{
    Real x,y,z,w;

    Quaternion()
        : x(0), y(0), z(0), w(1)
    {
    }

    Quaternion(const Vec4r& q)
        : x(q[0]), y(q[1]), z(q[2]), w(q[3])
    {
    }

    Quaternion(Real x_, Real y_, Real z_, Real w_)
        : x(x_), y(y_), z(z_), w(w_)
    {
    }

    static Quaternion Identity()
    {
        return Quaternion(0, 0, 0, 1);
    }

    Quaternion conjugate() const
    {
        return Quaternion(-x, -y, -z, w);
    }

    Quaternion operator*(const Quaternion& q)
    {
        return Quaternion(
            w * q.x + x * q.w + y * q.z - z * q.y,  // x
            w * q.y - x * q.z + y * q.w + z * q.x,  // y
            w * q.z + x * q.y - y * q.x + z * q.w,  // z
            w * q.w - x * q.x - y * q.y - z * q.z   // w
        );
    }

    void normalize()
    {
        Real n = std::sqrt(x*x + y*y + z*z + w*w);
        x /= n;
        y /= n;
        z /= n;
        w /= n;
    }

    
};