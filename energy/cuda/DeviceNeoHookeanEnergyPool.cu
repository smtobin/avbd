#include "energy/cuda/DeviceNeoHookeanEnergyPool.cuh"
#include "energy/elastic/NeoHookeanEnergyPool.hpp"

namespace Energy
{

void DeviceNeoHookeanEnergyPool::upload(const NeoHookeanEnergyPool& src)
{
    DeviceTombstonePool::upload(src);
    data.upload(src.data);
}

} // namespace Energy