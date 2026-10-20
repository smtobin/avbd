#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "common/cuda/Vec6.cuh"
#include "energy/elastic/NeoHookeanEnergyPool.hpp"

namespace Energy
{

struct DeviceCosseratRodEnergyInfo
{
    Vec2u particle_indices;
    DVec6 stiffness;
    float rest_length;
    DVec3 precurvature; 
};

struct DeviceCosseratRodEnergyPoolView
{
    DeviceCosseratRodEnergyInfo* data;
};

struct DeviceCosseratRodEnergyPool
{
    DeviceBuffer<DeviceCosseratRodEnergyInfo> data;
    
    explicit DeviceCosseratRodEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {

    }

    DeviceCosseratRodEnergyPoolView view() const
    {
        DeviceCosseratRodEnergyPoolView v;
        DeviceTombstonePoolView::view(v);
        v.data = data.data();

        return v;
    }

    void upload(const CosseratRodEnergyPool& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data);
    }
};

} // namespace Energy