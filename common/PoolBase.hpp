#pragma once

#include "common/common.hpp"

#include <vector>

/** Base struct for all memory pools in the sim. */
struct PoolBase
{
    /** Total capacity of the pool */
    unsigned capacity;
    /** Number of active elements in the pool */
    unsigned size;

    PoolBase(unsigned capacity_)
        : capacity(capacity_)
        , size(0)
    {}
};