#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

template <typename EnergyInfo, int N>
struct DeviceHardConstraintEnergyPoolView : DeviceTombstonePoolView
{
    EnergyInfo* data;
};

template <typename EnergyInfo, int N>
struct DeviceHardConstraintEnergyPool : DeviceTombstonePool
{
    DeviceBuffer<EnergyInfo> data;

    explicit DeviceHardConstraintEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {}

    void view(DeviceHardConstraintEnergyPoolView<EnergyInfo, N>& v) const
    {
        DeviceTombstonePool::view(v);
        v.data = data.data();
    }

    void upload(const HardConstraintEnergyPool<EnergyInfo, N>& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};

} // namespace Energy