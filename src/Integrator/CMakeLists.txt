set(Integrator_inst
    CheapSynchronousRK2
    Integrator
   )


set(Integrator_sources )

instantiate(Integrator_inst Integrator_sources)

set(Integrator_headers
    CheapSynchronousRK2.hh
    Integrator.hh
    IntegratorInline.hh
    )

blt_add_library(NAME        Spheral_Integrator
                SOURCES     ${Integrator_sources}
                HEADERS     ${Integrator_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Integrator
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Integrator_headers}
        DESTINATION include/Integrator
        )
