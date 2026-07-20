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

install(TARGETS Common EXPORT Common-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${COMMON_HEADERS}
    DESTINATION include/${PROJECT_NAME}/common)