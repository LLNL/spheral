set(FieldOperations_inst
    sampleMultipleFields2Lattice
    sampleMultipleFields2LatticeMash
   )


set(FieldOperations_sources )

instantiate(FieldOperations_inst FieldOperations_sources)

set(FieldOperations_headers
    sampleMultipleFields2Lattice.hh
    )

blt_add_library(NAME        Spheral_FieldOperations
                SOURCES     ${FieldOperations_sources}
                HEADERS     ${FieldOperations_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_FieldOperations
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FieldOperations_headers}
        DESTINATION include/FieldOperations
        )
