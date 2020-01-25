set(NodeGenerators_inst
   )

set(NodeGenerators_sources
    fillFacetedVolume.cc
    generateCylDistributionFromRZ.cc
    )

instantiate(NodeGenerators_inst NodeGenerators_sources)

set(NodeGenerators_headers
    fillFacetedVolume.hh
    generateCylDistributionFromRZ.hh
    )

blt_add_library(NAME        Spheral_NodeGenerators
                SOURCES     ${NodeGenerators_sources}
                HEADERS     ${NodeGenerators_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_NodeGenerators
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${NodeGenerators_headers}
        DESTINATION include/NodeGenerators
        )
