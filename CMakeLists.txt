cmake_minimum_required(VERSION 3.13)
project(XPBDRods LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)                       # specify C++ standard
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-ftime-trace)
    message(STATUS "Added -ftime-trace flag")
endif()

find_package(Eigen3 3.3 REQUIRED)   # find Eigen
find_package(VTK COMPONENTS
    CommonColor
    CommonCore
    FiltersSources
    IOImage
    InteractionStyle
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Unable to find the VTK build folder.")
endif()

# find assimp
find_package(assimp REQUIRED)

include_directories("./")

# add_library(Common)
# add_subdirectory("common")
# target_link_libraries(Common Eigen3::Eigen assimp)

# add_library(Collision)
# add_subdirectory("collision")
# target_link_libraries(Collision Eigen3::Eigen Common)

# add_library(Solver)
# add_subdirectory("solver")
# target_link_libraries(Solver Eigen3::Eigen Common)

# add_library(Constraint)
# add_subdirectory("constraint")
# target_link_libraries(Constraint Eigen3::Eigen Common SimObject)

# add_library(SimObject)
# add_subdirectory("simobject")
# target_link_libraries(SimObject Eigen3::Eigen Common Constraint)

add_library(Graphics)
add_subdirectory("graphics")
target_link_libraries(Graphics Eigen3::Eigen ${VTK_LIBRARIES})

add_library(Simulation)
add_subdirectory("simulation")
target_link_libraries(Simulation Eigen3::Eigen ${VTK_LIBRARIES} Graphics yaml-cpp)


#####################################################
# SIMULATION EXECUTABLES
#####################################################
add_executable(Sim exec/main.cpp)
target_link_libraries(Sim Simulation ${VTK_LIBRARIES})





# vtk_module_autoinit is needed

vtk_module_autoinit(
  TARGETS Sim
  MODULES ${VTK_LIBRARIES}
)

