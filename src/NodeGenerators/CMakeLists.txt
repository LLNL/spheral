set(NodeGenerators_inst
  relaxNodeDistribution
  centroidalRelaxNodesImpl
  compactFacetedVolumes
  chooseRandomNonoverlappingCenter
  )

set(NodeGenerators_sources
  generateCylDistributionFromRZ.cc
  fillFacetedVolume.cc
  #readSiloPolyMesh.cc
  )
instantiate(NodeGenerators_inst NodeGenerators_sources)

set(NodeGenerators_headers
  centroidalRelaxNodesImpl.hh
  chooseRandomNonoverlappingCenter.hh
  compactFacetedVolumes.hh
  fillFacetedVolume.hh
  generateCylDistributionFromRZ.hh
  #readSiloPolyMesh.hh
  relaxNodeDistribution.hh
  )

set(NodeGenerators_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/NodeGenerators
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_NodeGenerators
                HEADERS     ${NodeGenerators_headers}
                SOURCES     ${NodeGenerators_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${NodeGenerators_includes}
                )

install(TARGETS             Spheral_NodeGenerators
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${NodeGenerators_headers}
        DESTINATION include/NodeGenerators
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_NodeGenerators)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
