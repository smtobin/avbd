cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(GRAPHICS_DIR "./")

set(GRAPHICS_HEADERS
    ${GRAPHICS_DIR}/GraphicsObject.hpp
    ${GRAPHICS_DIR}/MeshGraphicsObject.hpp
    ${GRAPHICS_DIR}/GraphicsScene.hpp
    ${GRAPHICS_DIR}/CustomVTKInteractorStyle.hpp
)

set(GRAPHICS_SOURCES
    ${GRAPHICS_DIR}/MeshGraphicsObject.cpp
    ${GRAPHICS_DIR}/GraphicsScene.cpp
    ${GRAPHICS_DIR}/CustomVTKInteractorStyle.cpp
)

target_sources(Graphics
    PRIVATE
        ${GRAPHICS_HEADERS}
        ${GRAPHICS_SOURCES}
        
)

install(TARGETS Graphics EXPORT Graphics-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${GRAPHICS_HEADERS}
    DESTINATION include/${PROJECT_NAME}/graphics)