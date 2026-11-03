#pragma once

#include "common/common.hpp"
#include "graphics/VisualizationState.hpp"

namespace Graphics
{

struct VisualizationSlot
{
    VisualizationState host;
    cudaEvent_t copy_done;
    bool in_flight = false;
};

} // namespace Graphics