#pragma once

#include "energy/cuda/DeviceCollisionConstraintEnergyPool.cuh"
#include "energy/collision/GroundCollisionEnergyPool.hpp"

namespace Energy
{

struct DeviceGroundCollisionEnergyPoolView : DeviceCollisionConstraintEnergyPoolView<GroundCollisionEnergyInfo>
{
};

struct DeviceGroundCollisionEnergyPool : DeviceCollisionConstraintEnergyPool<GroundCollisionEnergyInfo>
{
    explicit DeviceGroundCollisionEnergyPool(unsigned capacity)
        :  DeviceCollisionConstraintEnergyPool(capacity)
    {
    }

    DeviceGroundCollisionEnergyPoolView view() const
    {
        DeviceGroundCollisionEnergyPoolView v;
        DeviceCollisionConstraintEnergyPool::view(v);
        return v;
    }

    void upload(const GroundCollisionEnergyPool& src)
    {
        DeviceCollisionConstraintEnergyPool::upload(src);
    }
};

} // namespace Energy