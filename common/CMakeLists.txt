cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(COMMON_DIR "./")

set(COMMON_HEADERS
    ${COMMON_DIR}/mesh/ParticleMesh.hpp
    ${COMMON_DIR}/mesh/ParticleTetMesh.hpp
    ${COMMON_DIR}/mesh/MeshUtils.hpp
)

set(COMMON_SOURCES
    ${COMMON_DIR}/mesh/ParticleMesh.cpp
    ${COMMON_DIR}/mesh/ParticleTetMesh.cpp
    ${COMMON_DIR}/mesh/MeshUtils.cpp
)

target_sources(Common
    PRIVATE
        ${COMMON_HEADERS}
        ${COMMON_SOURCES}
        
)

if (ENABLE_CUDA)
    # Example
    # set(SIMULATION_CUDA_SOURCES
    #     ${SIMULATION_DIR}/cuda/CudaSolver.cu
    # )

    # target_sources(Simulation
    #     PRIVATE
    #         ${SIMULATION_CUDA_SOURCES}
    # )

    target_link_libraries(Common CUDA::cudart)
    target_compile_definitions(Common PUBLIC HAS_CUDA=1)
else ()
    target_compile_definitions(Common PUBLIC HAS_CUDA=0)
endif()

install(TARGETS Common EXPORT Common-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${COMMON_HEADERS}
    DESTINATION include/${PROJECT_NAME}/common)