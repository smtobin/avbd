#pragma once

#include "energy/cuda/DeviceHardConstraintEnergyPool.cuh"
#include "energy/joint/OneSidedFixedJointEnergyPool.hpp"

namespace Energy
{

struct DeviceOneSidedFixedJointEnergyPool : DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>
{
    explicit DeviceOneSidedFixedJointEnergyPool(unsigned capacity)
        : DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>(capacity)
    {

    } 

    OneSidedFixedJointEnergyPoolView view()
    {
        OneSidedFixedJointEnergyPoolView v;
        DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>::view(v);
        return v;
    }

    void upload(const OneSidedFixedJointEnergyPool& src)
    {
        DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>::upload(src);
    }
};

} // namespace Energy