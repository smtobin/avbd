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

    static DeviceNeoHookeanEnergyInfo from(const NeoHookeanEnergyInfo& src)
    {
        DeviceNeoHookeanEnergyInfo i;
        i.particle_indices = src.particle_indices;
        i.lambda = src.lambda;
        i.mu = src.mu;
        i.kd = src.kd;
        i.Q = DMat3::from(Q);
        i.rest_volume = src.rest_volume;
        i.E_prev = DMat3::from(E_prev);
    }
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
        DeviceTombstonePool::upload(src);
        data.upload(src.data);
    }
};

} // namespace Energy