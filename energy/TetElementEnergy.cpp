#include "energy/TetElementEnergy.hpp"

namespace Energy
{

TetElementEnergy::TetElementEnergy(const std::array<const Particle*, 4>& element_particles)
    : _particles(element_particles)
{

}

Real TetElementEnergy::energy() const
{
    
}

} // namespace Energy