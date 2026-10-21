#pragma once

#include "energy/cuda/DeviceHardConstraintEnergyPool.cuh"
#include "energy/joint/OneSidedFixedJointEnergyPool.hpp"

namespace Energy
{

struct DeviceOneSidedFixedJointEnergyPoolView : DeviceHardConstraintEnergyPoolView<OneSidedFixedJointEnergyInfo, 6>
{
};

struct DeviceOneSidedFixedJointEnergyPool : DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>
{
    explicit DeviceOneSidedFixedJointEnergyPool(unsigned capacity)
        : DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo, 6>(capacity)
    {

    } 

    DeviceOneSidedFixedJointEnergyPoolView view() const
    {
        DeviceOneSidedFixedJointEnergyPoolView v;
        DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo>::view(v);
        return v;
    }

    void upload(const OneSidedFixedJointEnergyPool& src)
    {
        DeviceHardConstraintEnergyPool<OneSidedFixedJointEnergyInfo>::upload(src);
    }
};

} // namespace Energy