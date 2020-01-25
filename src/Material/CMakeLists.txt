set(Material_inst
    EquationOfState
   )


set(Material_sources
    PhysicalConstants.cc
    )

instantiate(Material_inst Material_sources)

set(Material_headers
    EquationOfState.hh
    EquationOfStateInline.hh
    PhysicalConstants.hh
    PhysicalConstantsInline.hh
    )

blt_add_library(NAME        Spheral_Material
                SOURCES     ${Material_sources}
                HEADERS     ${Material_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Material
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Material_headers}
        DESTINATION include/Material
        )
