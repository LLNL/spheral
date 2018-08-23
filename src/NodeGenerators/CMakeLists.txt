set(NodeGenerators_inst
    centroidalRelaxNodesImpl
    chooseRandomNonoverlappingCenter
    compactFacetedVolumes
    relaxNodeDistribution
   )

set(NodeGenerators_sources
    generateCylDistributionFromRZ.cc
    fillFacetedVolume.cc
    readSiloPolyMesh.cc
    )

instantiate(NodeGenerators_inst NodeGenerators_sources)

set(NodeGenerators_headers
    centroidalRelaxNodesImpl.hh
    chooseRandomNonoverlappingCenter.hh
    compactFacetedVolumes.hh
    fillFacetedVolume.hh
    generateCylDistributionFromRZ.hh
    readSiloPolyMesh.hh
    relaxNodeDistribution.hh
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
