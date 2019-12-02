set(Gravity_inst
  CompatibleGravitationalVelocityPolicy
  NBodyGravity
  )

set(Gravity_sources
  #LaplacianMatrixFactory.cc
  TreeGravity.cc
  )
instantiate(Gravity_inst Gravity_sources)

set(Gravity_headers
  CompatibleGravitationalVelocityPolicy.hh
  #LaplacianMatrixFactory.hh
  NBodyGravity.hh
  TreeGravity.hh
  TreeGravityInline.hh
  )

set(Gravity_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Gravity
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Gravity
                HEADERS     ${Gravity_headers}
                SOURCES     ${Gravity_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Gravity_includes}
                SHARED      TRUE
                )

install(TARGETS             Spheral_Gravity
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Gravity_headers}
        DESTINATION include/Gravity
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Gravity)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
