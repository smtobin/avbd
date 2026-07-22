#pragma once

#include "common/common.hpp"

namespace Sim
{

/** Simple struct that stores simulation parameters. */
struct SimulationParams
{
    Real dt;
    Real end_time;
    Real g_accel;
    Real viewer_refresh_time_ms;
    Real collision_margin;
};

} // namespace Sim