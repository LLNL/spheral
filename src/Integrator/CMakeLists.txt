set(Integrator_inst
  Integrator
  PredictorCorrector
  SynchronousRK1
  SynchronousRK2
  SynchronousRK4
  CheapSynchronousRK2
  Verlet
  )

set(Integrator_sources )
instantiate(Integrator_inst Integrator_sources)

set(Integrator_headers
  CheapSynchronousRK2.hh
  Integrator.hh
  IntegratorInline.hh
  PredictorCorrector.hh
  SynchronousRK1.hh
  SynchronousRK2.hh
  SynchronousRK4.hh
  Verlet.hh
  )

set(Integrator_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Integrator
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Integrator
                HEADERS     ${Integrator_headers}
                SOURCES     ${Integrator_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Integrator_includes}
                )

install(TARGETS             Spheral_Integrator
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Integrator_headers}
        DESTINATION include/Integrator
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Integrator)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
