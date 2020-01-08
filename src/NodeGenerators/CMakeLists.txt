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

set(NodeGenerators_python_targets
    NodeGeneratorBase.py
    GenerateNodeDistribution2d.py
    GenerateNodeDistribution3d.py
    CompositeNodeDistribution.py
    distributeNodesGeneric.py
    DistributeNodes.py
    ParMETISDistributeNodes.py
    MortonOrderDistributeNodes.py
    PeanoHilbertDistributeNodes.py
    VoronoiDistributeNodes.py
    NestedGridDistributeNodes.py
    GzipFileNodeGenerator.py
    DumpGzipFileNodeGenerator.py
    CubicNodeGenerator.py
    GenerateSphericalShellSection.py
    AbaqusNodeGenerator.py
    CentroidalVoronoiRelaxation.py
    InteriorGenerator.py
    FacetedSurfaceGenerator.py
    AsciiFileNodeGenerator.py
    HydroStaticProfile.py
    StretchedLattice2d.py
    StretchedLattice3d.py
    GenerateEqualMassSheets3d.py
    SiloPolyMeshGenerator.py
    GenerateRatioSlab2d.py
    GenerateRatioSphere.py
    GenerateNodeProfile.py
    CentroidalGenerator2d.py
    centroidalRelaxNodes.py
    MedialGenerator.py
    MultiScaleMedialGenerator.py
    RecursivePrimitiveRefinementGenerator.py
    SEAGenerator.py
    )

spheral_add_cxx_library(NodeGenerators)

install(TARGETS              Spheral_NodeGenerators
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${NodeGenerators_headers}
        DESTINATION include/NodeGenerators
        )

if(NOT ENABLE_CXXONLY)
  install(FILES       ${NodeGenerators_python_targets}
          DESTINATION .
          )
endif()
list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_NodeGenerators)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
