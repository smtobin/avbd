#pragma once

#include "energy/OneSidedFixedJointEnergyPool.hpp"

namespace Energy
{

struct OneSidedFixedJointConstraintSolver
{
    using PoolType = OneSidedFixedJointEnergyPool;

    static Vec6r evaluateConstraint(
        unsigned c_idx,
        const OneSidedFixedJointEnergyPool& energies,
        ParticlePool& particles,
    )
    {
        /** TODO */
    }
}

    
} // namespace Energy