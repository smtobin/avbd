#pragma once

#include "common/cuda/Vec3.cuh"
#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/collision/CollisionConstraintEnergyPool.hpp"

namespace Energy
{

struct DeviceCollisionConstraintEnergyInfo
{
    float k_start;
    float k_n; float k_t; float k_b;
    float lambda_n; float lambda_t; float lambda_b;
    float C_n_prev;
    DVec3 normal;
    DVec3 tangent;
    DVec3 binormal;

    bool use_static;
    bool unclamped_lambda_tb;
    float mu_s;
    float mu_k;

    static DeviceNeoHoo
};



} // namespace Energy