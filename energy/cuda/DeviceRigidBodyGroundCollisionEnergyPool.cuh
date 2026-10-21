#pragma once

#include "energy/cuda/DeviceCollisionConstraintEnergyPool.cuh"
#include "energy/collision/RigidBodyGroundCollisionEnergyPool.hpp"

namespace Energy
{

struct DeviceRigidBodyGroundCollisionEnergyPoolView : DeviceCollisionConstraintEnergyPoolView<RigidBodyGroundCollisionEnergyInfo>
{
};

struct DeviceRigidBodyGroundCollisionEnergyPool : DeviceCollisionConstraintEnergyPool<RigidBodyGroundCollisionEnergyInfo>
{
    explicit DeviceRigidBodyGroundCollisionEnergyPool(unsigned capacity)
        :  DeviceCollisionConstraintEnergyPool(capacity)
    {
    }

    DeviceRigidBodyGroundCollisionEnergyPoolView view() const
    {
        DeviceRigidBodyGroundCollisionEnergyPoolView v;
        DeviceCollisionConstraintEnergyPool::view(v);
        return v;
    }

    void upload(const RigidBodyGroundCollisionEnergyPool& src)
    {
        DeviceCollisionConstraintEnergyPool::upload(src);
    }
};

} // namespace Energy