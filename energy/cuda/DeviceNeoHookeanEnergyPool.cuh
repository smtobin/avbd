#pragma once

#include "common/cuda/Mat3.cuh"
#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/elastic/NeoHookeanEnergyPool.hpp"

namespace Energy
{

struct DeviceNeoHookeanEnergyInfo
{
    Vec4u particle_indices;
    float lambda;
    float mu;
    float kd;
    DMat3 Q;
    float rest_volume;
    DMat3 E_prev;

    /** TODO: (10/19/26) For custom device types, have a uniform method to convert from analogous CPU type. */
};

struct DeviceNeoHookeanEnergyPoolView : DeviceTombstonePoolView
{
    DeviceNeoHookeanEnergyInfo* data;
};

struct DeviceNeoHookeanEnergyPool : DeviceTombstonePool
{
    DeviceBuffer<NeoHookeanEnergyInfo> data;
    
    explicit DeviceNeoHookeanEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {

    }

    DeviceNeoHookeanEnergyPoolView view() const
    {
        DeviceNeoHookeanEnergyPoolView v;
        DeviceTombstonePoolView::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const NeoHookeanEnergyPool& src)
    {
        /** TODO: (10/19/26) Finish this. */
    }
};

} // namespace Energy