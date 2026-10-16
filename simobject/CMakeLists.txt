cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(SIMOBJECT_DIR "./")

set(SIMOBJECT_HEADERS
    ${SIMOBJECT_DIR}/ObjectBase.hpp
    ${SIMOBJECT_DIR}/TetMeshObject.hpp
    ${SIMOBJECT_DIR}/Rod.hpp

    ${SIMOBJECT_DIR}/rigid/RigidObjectBase.hpp
    ${SIMOBJECT_DIR}/rigid/RigidSphere.hpp
)

set(SIMOBJECT_SOURCES
    ${SIMOBJECT_DIR}/TetMeshObject.cpp
    ${SIMOBJECT_DIR}/Rod.cpp

    ${SIMOBJECT_DIR}/rigid/RigidObjectBase.cpp
    ${SIMOBJECT_DIR}/rigid/RigidSphere.cpp
)

target_sources(SimObject
    PRIVATE
        ${SIMOBJECT_HEADERS}
        ${SIMOBJECT_SOURCES}
        
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

    target_link_libraries(SimObject CUDA::cudart)
    target_compile_definitions(SimObject PUBLIC HAS_CUDA=1)
else ()
    target_compile_definitions(SimObject PUBLIC HAS_CUDA=0)
endif()

install(TARGETS SimObject EXPORT SimObject-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${SIMOBJECT_HEADERS}
    DESTINATION include/${PROJECT_NAME}/simobject)