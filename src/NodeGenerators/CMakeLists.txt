include_directories(.)
set(NodeGenerators_inst
    relaxNodeDistribution
    centroidalRelaxNodesImpl
    compactFacetedVolumes
    chooseRandomNonoverlappingCenter
   )

set(NodeGenerators_sources
    generateCylDistributionFromRZ.cc
    fillFacetedVolume.cc
    readSiloPolyMesh.cc
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
                SHARED      TRUE
                )

install(TARGETS              Spheral_NodeGenerators
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${NodeGenerators_headers}
        DESTINATION include/NodeGenerators
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_NodeGenerators)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
