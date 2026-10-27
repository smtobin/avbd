#pragma once

#include "common/cuda/DeviceTombstonePool.cuh"
#include "energy/collision/CollisionConstraintEnergyDefs.hpp"
#include "energy/collision/CollisionConstraintEnergyPool.hpp"

namespace Energy
{

template <typename CollisionEnergyInfo>
struct DeviceCollisionConstraintEnergyPool
{
    DeviceBuffer<CollisionEnergyInfo> data;

    explicit DeviceCollisionConstraintEnergyPool(unsigned capacity)
        : DeviceTombstonePool(capacity)
        , data(capacity)
    {}

    void view(CollisionConstraintEnergyPoolView<CollisionEnergyInfo>& v) const
    {
        DeviceTombstonePoolView::view(v);
        v.data = data.data();
    }

    void upload(const CollisionConstraintEnergyPool<CollisionEnergyInfo>& src)
    {
        DeviceTombstonePool::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};



} // namespace Energy