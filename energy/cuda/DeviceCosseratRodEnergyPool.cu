#include "energy/cuda/DeviceCosseratRodEnergyPool.cuh"
#include "energy/elastic/CosseratRodEnergyPool.hpp"

namespace Energy
{

void DeviceCosseratRodEnergyPool::upload(const CosseratRodEnergyPool& src);
{
    DeviceTombstonePool::upload(src);
    data.upload(src.data);
}

} // namespace Energy