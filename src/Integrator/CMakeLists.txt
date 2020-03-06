include_directories(.)
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
    )

spheral_add_cxx_library(Integrator)

install(TARGETS              Spheral_Integrator
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Integrator_headers}
        DESTINATION include/Integrator
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Integrator)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
