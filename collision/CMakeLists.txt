cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(COLLISION_DIR "./")

set(COLLISION_HEADERS
    ${COLLISION_DIR}/CollisionObject.hpp
    ${COLLISION_DIR}/CollisionDetector.hpp
    ${COLLISION_DIR}/CollisionPrimitivePool.hpp
    ${COLLISION_DIR}/LBVH.hpp
    ${COLLISION_DIR}/LBVHBuilder.hpp
    ${COLLISION_DIR}/LBVHTraversal.hpp
    ${COLLISION_DIR}/AABB.hpp
)

set(COLLISION_SOURCES
    ${COLLISION_DIR}/CollisionObject.cpp
    ${COLLISION_DIR}/CollisionDetector.cpp
    ${COLLISION_DIR}/CollisionPrimitivePool.cpp
    ${COLLISION_DIR}/LBVH.cpp
    ${COLLISION_DIR}/LBVHBuilder.cpp
)

target_sources(Collision
    PRIVATE
        ${COLLISION_HEADERS}
        ${COLLISION_SOURCES}
        
)

install(TARGETS Collision EXPORT Collision-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${COLLISION_HEADERS}
    DESTINATION include/${PROJECT_NAME}/simobject)