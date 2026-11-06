#pragma once

#include "common/cuda/DevicePoolBase.cuh"
#include "energy/HardConstraintEnergyPool.hpp"

namespace Energy
{

template <typename EnergyInfo, int N>
struct DeviceHardConstraintEnergyPool : public DevicePoolBase
{
    DeviceBuffer<EnergyInfo> data;

    explicit DeviceHardConstraintEnergyPool(unsigned capacity)
        : DevicePoolBase(capacity)
        , data(capacity)
    {}

    void view(HardConstraintEnergyPoolView<EnergyInfo, N>& v)
    {
        DevicePoolBase::view(v);
        v.data = data.data();
    }

    void upload(const HardConstraintEnergyPool<EnergyInfo, N>& src)
    {
        DevicePoolBase::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};

} // namespace Energy