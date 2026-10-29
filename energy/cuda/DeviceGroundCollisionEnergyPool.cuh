#pragma once

#include "energy/cuda/DeviceCollisionConstraintEnergyPool.cuh"
#include "energy/collision/GroundCollisionEnergyPool.hpp"

namespace Energy
{

struct DeviceGroundCollisionEnergyPool : DeviceCollisionConstraintEnergyPool<GroundCollisionEnergyInfo>
{
    explicit DeviceGroundCollisionEnergyPool(unsigned capacity)
        :  DeviceCollisionConstraintEnergyPool<GroundCollisionEnergyInfo>(capacity)
    {
    }

    GroundCollisionEnergyPoolView view()
    {
        GroundCollisionEnergyPoolView v;
        DeviceCollisionConstraintEnergyPool<GroundCollisionEnergyInfo>::view(v);
        return v;
    }

    void upload(const GroundCollisionEnergyPool& src)
    {
        DeviceCollisionConstraintEnergyPool<GroundCollisionEnergyInfo>::upload(src);
    }
};

} // namespace Energy