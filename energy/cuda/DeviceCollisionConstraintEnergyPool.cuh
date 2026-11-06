#pragma once

#include "common/cuda/DevicePoolBase.cuh"
#include "energy/collision/CollisionConstraintEnergyDefs.hpp"
#include "energy/collision/CollisionConstraintEnergyPool.hpp"

namespace Energy
{

template <typename CollisionEnergyInfo>
struct DeviceCollisionConstraintEnergyPool : public DevicePoolBase
{
    DeviceBuffer<CollisionEnergyInfo> data;

    explicit DeviceCollisionConstraintEnergyPool(unsigned capacity)
        : DevicePoolBase(capacity)
        , data(capacity)
    {}

    void view(CollisionConstraintEnergyPoolView<CollisionEnergyInfo>& v)
    {
        DevicePoolBase::view(v);
        v.data = data.data();
    }

    void upload(const CollisionConstraintEnergyPool<CollisionEnergyInfo>& src)
    {
        DevicePoolBase::upload(src);
        data.upload(src.data.data(), src.data.size());
    }
};



} // namespace Energy