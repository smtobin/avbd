cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(COMMON_DIR "./")

set(COMMON_HEADERS
    ${COMMON_DIR}/mesh/Mesh.hpp
    ${COMMON_DIR}/mesh/TetMesh.hpp
)

set(COMMON_SOURCES
    ${COMMON_DIR}/mesh/Mesh.cpp
    ${COMMON_DIR}/mesh/TetMesh.cpp
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
    DESTINATION include/${PROJECT_NAME}/graphics)