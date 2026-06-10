cmake_minimum_required(VERSION 3.13)
project(XPBDRods)

set(ENERGY_DIR "./")

set(ENERGY_HEADERS
    # energies
    ${ENERGY_DIR}/EnergyBase.hpp
    ${ENERGY_DIR}/QuadraticEnergy.hpp
    ${ENERGY_DIR}/TetElementEnergy.hpp
    ${ENERGY_DIR}/GroundCollisionEnergy.hpp
    ${ENERGY_DIR}/HardConstraintEnergy.hpp

    # constraints
    ${ENERGY_DIR}/constraint/ConstraintBase.hpp
    ${ENERGY_DIR}/constraint/AttachmentConstraint.hpp
    ${ENERGY_DIR}/constraint/GroundCollisionConstraint.hpp
)

set(ENERGY_SOURCES
    # energies
    ${ENERGY_DIR}/EnergyBase.cpp
    ${ENERGY_DIR}/TetElementEnergy.cpp
    ${ENERGY_DIR}/GroundCollisionEnergy.cpp
    ${ENERGY_DIR}/HardConstraintEnergy.cpp

    # constraints
    ${ENERGY_DIR}/constraint/AttachmentConstraint.cpp
    ${ENERGY_DIR}/constraint/GroundCollisionConstraint.cpp
)

target_sources(Energy
    PRIVATE
        ${ENERGY_HEADERS}
        ${ENERGY_SOURCES}
        
)

install(TARGETS Energy EXPORT Energy-targets
    RUNTIME DESTINATION bin/${PROJECT_NAME}
    LIBRARY DESTINATION lib/${PROJECT_NAME})
install(FILES ${ENERGY_HEADERS}
    DESTINATION include/${PROJECT_NAME}/graphics)