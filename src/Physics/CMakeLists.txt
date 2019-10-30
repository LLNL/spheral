set(Physics_inst
  GenericBodyForce
  GenericHydro
  Physics)

set(Physics_sources )
instantiate(Physics_inst Physics_sources)

set(Physics_headers
  GenericBodyForce.hh
  GenericHydro.hh
  Physics.hh
  )

set(Physics_includes
  ${SPHERAL_SOURCE_DIR}
  ${EIGEN_LIB_LOCATION}
  )

blt_add_library(NAME Spheral_Physics
  HEADERS    ${Physics_headers}
  SOURCES    ${Physics_sources}
  DEPENDS_ON ${spheral_depends}
  INCLUDES   ${Physics_includes}
  )

install(TARGETS Spheral_Physics
  EXPORT               ${PROJECT_NAME}-targets
  DESTINATION          lib
  INCLUDES DESTINATION include
  )

install(FILES ${Physics_headers}
  DESTINATION include/Physics
  )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Physics)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
