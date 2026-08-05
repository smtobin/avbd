#pragma once

#include "common/common.hpp"

namespace Energy
{

/** Helper for iterating over different energies.
 */
// Usage:
//    ForEachEnergy([&]<EnergyType E>() {
//              // some work
//      });
template<std::size_t... Is, typename F>
constexpr void ForEachEnergy_Impl(std::index_sequence<Is...>, F&& f)
{
    (std::forward<F>(f).template operator()<static_cast<EnergyType>(Is)>(), ...);
}

template<typename F>
constexpr void ForEachEnergy(F&& f)
{
    ForEachEnergy_Impl(
        std::make_index_sequence<static_cast<std::size_t>(EnergyType::count)>{}, 
        std::forward<F>(f)    
    );
}

/** Helper for iterating over different static energies.
 */
// Usage:
//    ForEachStaticEnergy([&]<StaticEnergyType E>() {
//              // some work
//      });
template<std::size_t... Is, typename F>
constexpr void ForEachStaticEnergy_Impl(std::index_sequence<Is...>, F&& f)
{
    (std::forward<F>(f).template operator()<static_cast<StaticEnergyType>(Is)>(), ...);
}

template<typename F>
constexpr void ForEachStaticEnergy(F&& f)
{
    ForEachStaticEnergy_Impl(
        std::make_index_sequence<static_cast<std::size_t>(StaticEnergyType::count)>{}, 
        std::forward<F>(f)    
    );
}

/** Helper for iterating over different dynamic energies.
 */
// Usage:
//    ForEachStaticEnergy([&]<StaticEnergyType E>() {
//              // some work
//      });
template<std::size_t... Is, typename F>
constexpr void ForEachDynamicEnergy_Impl(std::index_sequence<Is...>, F&& f)
{
    (std::forward<F>(f).template operator()<static_cast<DynamicEnergyType>(Is)>(), ...);
}

template<typename F>
constexpr void ForEachDynamicEnergy(F&& f)
{
    ForEachDynamicEnergy_Impl(
        std::make_index_sequence<static_cast<std::size_t>(DynamicEnergyType::count)>{}, 
        std::forward<F>(f)    
    );
}



/** EnergySolver concepts */
template<typename Solver>
concept HasAccumulate4 = requires(
    unsigned e_idx[4],
    const Solver::PoolType& energies,
    ParticlePool& particles,
    unsigned local_idx[4],
    Mat3r& particle_H,
    Vec3r& particle_G,
    Real dt
)
{
    Solver::accumulate4(e_idx, energies, particles, local_idx, particle_H, particle_G, dt);
};

}