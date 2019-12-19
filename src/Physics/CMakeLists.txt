set(Physics_inst
    GenericBodyForce
    GenericHydro
    Physics
   )

set(Physics_sources )

instantiate(Physics_inst Physics_sources)

set(Physics_headers
    GenericBodyForce.hh
    GenericHydro.hh
    GenericHydroInline.hh
    Physics.hh
    PhysicsInline.hh
    )

blt_add_library(NAME        Spheral_Physics
                SOURCES     ${Physics_sources}
                HEADERS     ${Physics_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Physics
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Physics_headers}
        DESTINATION include/Physics
        )
