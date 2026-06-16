#pragma once

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

using Vec3u = Eigen::Vector<unsigned, 3>;
using Vec4u = Eigen::Vector<unsigned, 4>;

using Mat2r = Eigen::Matrix<Real, 2, 2>;
using Mat3r = Eigen::Matrix<Real, 3, 3>;
using Mat4r = Eigen::Matrix<Real, 4, 4>;
using Mat6r = Eigen::Matrix<Real, 6, 6>;
using MatXr = Eigen::Matrix<Real,-1,-1>;

/** Enum of energy types */
enum class EnergyType
{
    NEO_HOOKEAN = 0,
    GROUND_COLLISION
};

/** Forward declarations */
namespace Energy
{
    class Energy_Base;
}

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

using ObjectConfigs_TypeList = TypeList<
    Config::RigidSphereConfig,
    Config::TetMeshObjectConfig
>;

using ObjectConfigs_Container = VariadicVectorContainerFromTypeList<ObjectConfigs_TypeList>::type;