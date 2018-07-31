set(NodeGenerators_sources
    generateCylDistributionFromRZ.cc
    fillFacetedVolume.cc
    readSiloPolyMesh.cc
    )

set(NodeGenerators_sources2d
    centroidalRelaxNodesImplInst.cc.2.cc
    chooseRandomNonoverlappingCenterInst.cc.2.cc
    compactFacetedVolumesInst.cc.2.cc
    relaxNodeDistributionInst.cc.2.cc
    )

set(NodeGenerators_sources3d
    centroidalRelaxNodesImplInst.cc.3.cc
    chooseRandomNonoverlappingCenterInst.cc.3.cc
    compactFacetedVolumesInst.cc.3.cc
    relaxNodeDistributionInst.cc.3.cc
    )

set(NodeGenerators_headers
    centroidalRelaxNodesImpl.hh
    chooseRandomNonoverlappingCenter.hh
    compactFacetedVolumes.hh
    fillFacetedVolume.hh
    generateCylDistributionFromRZ.hh
    readSiloPolyMesh.hh
    relaxNodeDistribution.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_NodeGenerators
                SOURCES     ${NodeGenerators_sources} ${NodeGenerators_sources2d} ${NodeGenerators_sources3d}
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
