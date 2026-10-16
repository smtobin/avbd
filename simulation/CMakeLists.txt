cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(SIMULATION_DIR "./")

set(SIMULATION_HEADERS
    ${SIMULATION_DIR}/Simulation.hpp
    # ${SIMULATION_DIR}/SimulationLogger.hpp
)

set(SIMULATION_SOURCES
    ${SIMULATION_DIR}/Simulation.cpp
    # ${SIMULATION_DIR}/SimulationLogger.cpp
)

target_sources(Simulation
    PRIVATE
        ${SIMULATION_HEADERS}
        ${SIMULATION_SOURCES}
        
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

    target_link_libraries(Simulation CUDA::cudart)
    target_compile_definitions(Simulation PUBLIC HAS_CUDA=1)
else ()
    target_compile_definitions(Simulation PUBLIC HAS_CUDA=0)
endif()

install(TARGETS Simulation EXPORT Simulation-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${SIMULATION_HEADERS}
    DESTINATION include/${PROJECT_NAME}/simulation)