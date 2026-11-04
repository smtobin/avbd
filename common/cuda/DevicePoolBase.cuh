#pragma once

#include "common/common.hpp"
#include "common/PoolBase.hpp"
#include "common/cuda/DeviceValue.cuh"
#include "common/cuda/DeviceBuffer.cuh"

/** Append-only base for all GPU memory pools in the sim. */
struct DevicePoolBase
{
    unsigned capacity;
    unsigned host_high_water;

    DeviceValue<unsigned> high_water;
};