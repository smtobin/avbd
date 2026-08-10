#pragma once

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <memory>
#include <cassert>

#include "common/TypeList.hpp"
// #include "common/VariadicVectorContainer.hpp"

#define STIFFNESS_BETA 10
#define STIFFNESS_GAMMA 0.99
#define CONSTRAINT_ALPHA 0.95

/** Forward declaration of VariadicVectorContainer */
template<class L, class... R> class VariadicVectorContainer;

//////////////////////////////////////////////////////////////////////////
// Construct VariadicVectorContainer from TypeList
//////////////////////////////////////////////////////////////////////////

template<typename List>
struct VariadicVectorContainerFromTypeList;

template<typename... Types>
struct VariadicVectorContainerFromTypeList<TypeList<Types...>>
{
    using type = VariadicVectorContainer<Types...>;
    using unique_ptr_type = VariadicVectorContainer<std::unique_ptr<Types>...>;
    using ptr_type = VariadicVectorContainer<Types*...>;
    using const_ptr_type = VariadicVectorContainer<const Types*...>;
    // using vector_handle_type = VariadicVectorContainer<VectorHandle<Types>...>;
    // using const_vector_handle_type = VariadicVectorContainer<ConstVectorHandle<Types>...>;
};

/////////////////////////////////////////////////////////////////////////
// Get base type (remove keywords, references, etc.)
/////////////////////////////////////////////////////////////////////////
template<typename T>
struct base_type { using type = T; };

template<typename T>
struct base_type<T*> : base_type<T> {};

template<typename T>
struct base_type<T&> : base_type<T> {};

template<typename T>
struct base_type<T&&> : base_type<T> {};

template<typename T>
struct base_type<const T> : base_type<T> {};

template<typename T>
struct base_type<volatile T> : base_type<T> {};

template<typename T>
using base_type_t = typename base_type<T>::type;


/** Escape sequences to set print colors */
#define RST  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define BOLD  "\x1B[1m"
#define UNDL  "\x1B[4m"


/** Universal typedefs used by the simulation */
using Real = double;

using Vec2r = Eigen::Vector<Real, 2>;
using Vec3r = Eigen::Vector<Real, 3>;
using Vec4r = Eigen::Vector<Real, 4>;
using Vec6r = Eigen::Vector<Real, 6>;
using VecXr = Eigen::Vector<Real, -1>;

using Vec3i = Eigen::Vector<int, 3>;
using Vec4i = Eigen::Vector<int, 4>;

using Vec1u = Eigen::Vector<unsigned, 1>;
using Vec2u = Eigen::Vector<unsigned, 2>;
using Vec3u = Eigen::Vector<unsigned, 3>;
using Vec4u = Eigen::Vector<unsigned, 4>;

using Mat2r = Eigen::Matrix<Real, 2, 2>;
using Mat3r = Eigen::Matrix<Real, 3, 3>;
using Mat4r = Eigen::Matrix<Real, 4, 4>;
using Mat6r = Eigen::Matrix<Real, 6, 6>;
using MatXr = Eigen::Matrix<Real,-1,-1>;

using Quaternion = Eigen::Quaternion<Real>;

// template for switching between Vec3r and Vec6r depending on some template parameter
template <int DOF>
using Vec3r_or_Vec6r = std::conditional_t<DOF == 6, Vec6r, Vec3r>;
template <int DOF>
using Mat3r_or_Mat6r = std::conditional_t<DOF == 6, Mat6r, Mat3r>;

/** Enum of energy types */
// enum class EnergyType
// {
//     NEO_HOOKEAN = 0,
//     GROUND_COLLISION,
//     size    // stores the number of different energies
// };

/** Forward declarations */
namespace Energy
{
    class Energy_Base;
    
    /** Pools */
    struct NeoHookeanEnergyPool;
    struct GroundCollisionEnergyPool;
    struct RigidBodyGroundCollisionEnergyPool;
    struct TriangleRigidCollisionEnergyPool;

    /** Solvers */
    struct NeoHookeanEnergySolver;
    struct GroundCollisionEnergySolver;
    struct RigidBodyGroundCollisionEnergySolver;
    struct TriangleRigidCollisionEnergySolver;
}


/** Using macros, define the EnergyType enum and a mapping to the corresponding Solver type. */
// when a new energy is added, we must update the list below
#define ENERGY_LIST(X) \
    X(NEO_HOOKEAN, NeoHookeanEnergySolver) \
    X(GROUND_COLLISION, GroundCollisionEnergySolver) \
    X(RIGID_BODY_GROUND_COLLISION, RigidBodyGroundCollisionEnergySolver) \
    X(TRIANGLE_RIGID_COLLISION, TriangleRigidCollisionEnergySolver)

// "Static" energies are those that are generally not added or removed throughout the course of the simulation
// (unless topology changes)
#define STATIC_ENERGY_LIST(X) \
    X(NEO_HOOKEAN, NeoHookeanEnergySolver) \
    X(GROUND_COLLISION, GroundCollisionEnergySolver) \
    X(RIGID_BODY_GROUND_COLLISION, RigidBodyGroundCollisionEnergySolver)

// "Dynamic" energies are those that are added often throughout the course of the simulation - e.g. most collision constraints
#define DYNAMIC_ENERGY_LIST(X) \
    X(TRIANGLE_RIGID_COLLISION, TriangleRigidCollisionEnergySolver)

// generate the enum for all energies
enum class EnergyType : uint8_t
{
#define X(name, solver) name,
    ENERGY_LIST(X)
#undef X
    count    // count = total number of energies
};

// generate the enum for static energies
enum class StaticEnergyType : uint8_t
{
#define X(name, solver) name,
    STATIC_ENERGY_LIST(X)
#undef X
    count   // count = number of static energies
};

// generate the enum for dynamic energies
enum class DynamicEnergyType : uint8_t
{
#define X(name, solver) name,
    DYNAMIC_ENERGY_LIST(X)
#undef X    
    count   // count = number of dynamic energies
};

// generate the mapping from energy type -> solver type
// Usage:
//   SolverFor<EnergyType::NEO_HOOKEAN>::type ==> Energy::NeoHookeanEnergySolver
//   SolverFor<StaticEnergyType::NEO_HOOKEAN>::type ==> Energy::NeoHookeanEnergySolver
//   SolverFor<DynamicEnergyType::TRIANGLE_RIGID_COLLISION>::type ==> TriangleRigidCollisionEnergySolver
template<auto>
struct SolverFor;

#define X(name, solver)            \
template<>                         \
struct SolverFor<EnergyType::name> \
{                                  \
    using type = Energy::solver;   \
};

ENERGY_LIST(X)
#undef X

// generate the mapping from static energy type -> solver type
#define X(name, solver)            \
template<>                         \
struct SolverFor<StaticEnergyType::name> \
{                                  \
    using type = Energy::solver;   \
};

STATIC_ENERGY_LIST(X)
#undef X

// generate the mapping from dynamic energy type -> solver type
#define X(name, solver)            \
template<>                         \
struct SolverFor<DynamicEnergyType::name> \
{                                  \
    using type = Energy::solver;   \
};

DYNAMIC_ENERGY_LIST(X)
#undef X


/** Forward declarations of objects in sim */
namespace SimObject
{
    class Object_Base;
    class TetMeshObject;
    class RigidSphere;
}

namespace Config
{
    class Config;
    class ObjectConfig;
    class TetMeshObjectConfig;
    class RigidSphereConfig;
}

namespace Collision
{
    struct CollisionPrimitivePool;
    struct LBVH;
}

namespace Sim
{
    struct SimulationContext;
}

using ObjectConfigs_TypeList = TypeList<
    Config::RigidSphereConfig,
    Config::TetMeshObjectConfig
>;

using ObjectConfigs_Container = VariadicVectorContainerFromTypeList<ObjectConfigs_TypeList>::type;