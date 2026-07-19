// #pragma once

// #include "common/common.hpp"

// struct Quaternion
// {
//     Real x,y,z,w;

//     Quaternion()
//         : x(0), y(0), z(0), w(1)
//     {
//     }

//     Quaternion(const Vec4r& q)
//         : x(q[0]), y(q[1]), z(q[2]), w(q[3])
//     {
//     }

//     Quaternion(Real x_, Real y_, Real z_, Real w_)
//         : x(x_), y(y_), z(z_), w(w_)
//     {
//     }

//     static Quaternion Identity()
//     {
//         return Quaternion(0, 0, 0, 1);
//     }

//     Quaternion conjugate() const
//     {
//         return Quaternion(-x, -y, -z, w);
//     }

//     inline Quaternion operator*(const Quaternion& q)
//     {
//         return Quaternion(
//             w * q.x + x * q.w + y * q.z - z * q.y,  // x
//             w * q.y - x * q.z + y * q.w + z * q.x,  // y
//             w * q.z + x * q.y - y * q.x + z * q.w,  // z
//             w * q.w - x * q.x - y * q.y - z * q.z   // w
//         );
//     }

//     inline Vec3r operator*(const Vec3r& v)
//     {
//         Vec3r qv(x, y, z);

//         Vec3r t = 2.0 * qv.cross(v);
//         return v + w * t + qv.cross(t);
//     }

//     inline void normalize()
//     {
//         Real n = std::sqrt(x*x + y*y + z*z + w*w);
//         x /= n;
//         y /= n;
//         z /= n;
//         w /= n;
//     }

//     inline Mat3r toRotationMatrix() const
//     {
//         const Real tx = 2 * x;
//         const Real ty = 2 * y;
//         const Real tz = 2 * z;

//         const Real twx = tx * w;
//         const Real twy = ty * w;
//         const Real twz = tz * w;

//         const Real txx = tx * x;
//         const Real txy = ty * x;
//         const Real txz = tz * x;

//         const Real tyy = ty * y;
//         const Real tyz = tz * y;
//         const Real tzz = tz * z;

//         Mat3r R;

//         R(0,0) = 1 - (tyy + tzz);
//         R(0,1) = txy - twz;
//         R(0,2) = txz + twy;

//         R(1,0) = txy + twz;
//         R(1,1) = 1 - (txx + tzz);
//         R(1,2) = tyz - twx;

//         R(2,0) = txz - twy;
//         R(2,1) = tyz + twx;
//         R(2,2) = 1 - (txx + tyy);

//         return R;
//     }

    
// };