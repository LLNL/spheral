set(Utilities_inst
    boundingBox
    computeShepardsInterpolation
    globalBoundingVolumes
    integrateThroughMeshAlongSegment
    iterateIdealH
    medianPosition
    mortonOrderIndices
    nodeBoundingBoxes
    nodeOrdering
    numberDensity
    peanoHilbertOrderIndices
   )

set(Utilities_sources
    DBC.cc
    Process.cc
    KeyTraits.cc
    PeanoHilbertTransform.cc
    RegisterMPIDataTypes.cc
    RedistributionRegistrar.cc
    vectorBoundingBoxInst.cc
    lineSegmentIntersections.cc
    segmentIntersectEdges.cc
    pointOnPolygon.cc
    pointInPolygon.cc
    pointOnPolyhedron.cc
    pointInPolyhedron.cc
    refinePolyhedron.cc
    coarsenBinnedValuesInst.cc
    SpheralTimers.cc
    SpheralTimersAlternate.cc
    Timer.cc
    )

instantiate(Utilities_inst Utilities_sources)

set(Utilities_headers
    Bessel.hh
    CounterClockwiseComparator.hh
    DBC.hh
    DataTypeTraits.hh
    FastMath.hh
    FieldDataTypeTraits.hh
    Functors.hh
    GeometricUtilities.hh
    KeyTraits.hh
    NodeListHashMap.hh
    PairComparisons.hh
    PeanoHilbertTransform.hh
    Process.hh
    RedistributionNotification.hh
    RedistributionNotificationHandle.hh
    RedistributionRegistrar.hh
    RedistributionRegistrarInline.hh
    RegisterMPIDataTypes.hh
    SafeIndexMap.hh
    SpheralFunctions.hh
    Timer.hh
    Tree.hh
    TreeInline.hh
    allReduce.hh
    bisectRoot.hh
    bisectSearch.hh
    boundPointWithinBox.hh
    boundingBox.hh
    coarsenBinnedValues.hh
    comparisons.hh
    computeConvexHull.hh
    computeShepardsInterpolation.hh
    erff.hh
    fakestl.hh
    globalBoundingVolumes.hh
    globalNodeIDs.hh
    globalNodeIDsInline.hh
    hashes.hh
    initializeTau.hh
    integrateThroughMeshAlongSegment.hh
    intpow2.hh
    iterateIdealH.hh
    lineSegmentIntersections.hh
    log2.hh
    medianPosition.hh
    mortonOrderIndices.hh
    mpiUtilities.hh
    newtonRaphson.hh
    nodeBoundingBoxes.hh
    nodeBoundingBoxesInline.hh
    nodeOrdering.hh
    numberDensity.hh
    overlayRemapFields.hh
    packElement.hh
    peanoHilbertOrderIndices.hh
    planarReflectingOperator.hh
    pointDistances.hh
    pointInPolygon.hh
    pointInPolyhedron.hh
    pointOnPolygon.hh
    pointOnPolyhedron.hh
    refinePolyhedron.hh
    registerWithRedistribution.hh
    removeElements.hh
    rotationMatrix.hh
    safeInv.hh
    segmentIntersectEdges.hh
    simpsonsIntegration.hh
    sort_permutation.hh
    testBoxIntersection.hh
    timingUtilities.hh
    )

INCLUDE_DIRECTORIES(. ../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Utilities
                SOURCES     ${Utilities_sources}
                HEADERS     ${Utilities_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Utilities
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Utilities_headers}
        DESTINATION include/Utilities
        )
