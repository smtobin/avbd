#pragma once

#include "energy/cuda/DeviceCollisionConstraintEnergyPool.cuh"
#include "energy/collision/RigidBodyGroundCollisionEnergyPool.hpp"

namespace Energy
{

struct DeviceRigidBodyGroundCollisionEnergyPool : DeviceCollisionConstraintEnergyPool<RigidBodyGroundCollisionEnergyInfo>
{
    explicit DeviceRigidBodyGroundCollisionEnergyPool(unsigned capacity)
        :  DeviceCollisionConstraintEnergyPool<RigidBodyGroundCollisionEnergyInfo>(capacity)
    {
    }

    RigidBodyGroundCollisionEnergyPoolView view()
    {
        RigidBodyGroundCollisionEnergyPoolView v;
        DeviceCollisionConstraintEnergyPool<RigidBodyGroundCollisionEnergyInfo>::view(v);
        return v;
    }

    void upload(const RigidBodyGroundCollisionEnergyPool& src)
    {
        DeviceCollisionConstraintEnergyPool<RigidBodyGroundCollisionEnergyInfo>::upload(src);
    }
};

} // namespace Energy