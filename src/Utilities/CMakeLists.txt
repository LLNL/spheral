<<<<<<< HEAD
include_directories(.)
set(Utilities_inst
    iterateIdealH
    mortonOrderIndices
    peanoHilbertOrderIndices
    nodeOrdering
    nodeBoundingBoxes
    boundingBox
    globalBoundingVolumes
    integrateThroughMeshAlongSegment
    numberDensity
    medianPosition
    computeShepardsInterpolation
    overlayRemapFields
=======
set(Utilities_inst
    boundingBox
    computeShepardsInterpolation
    globalBoundingVolumes
    iterateIdealH
    mortonOrderIndices
    nodeBoundingBoxes
    peanoHilbertOrderIndices
>>>>>>> origin/master
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
<<<<<<< HEAD
    refinePolyhedron.cc
    coarsenBinnedValuesInst.cc
    clipFacetedVolume.cc
    SpheralTimers.cc
    #SpheralTimersAlternate.cc
=======
    coarsenBinnedValuesInst.cc
    SpheralTimers.cc
    SpheralTimersAlternate.cc
>>>>>>> origin/master
    Timer.cc
    )

instantiate(Utilities_inst Utilities_sources)

set(Utilities_headers
    Bessel.hh
    CounterClockwiseComparator.hh
    DBC.hh
<<<<<<< HEAD
    DomainNode.hh
    DomainNodeInline.hh
    DataTypeTraits.hh
=======
    DataTypeTraits.hh
    DomainNode.hh
    DomainNodeInline.hh
>>>>>>> origin/master
    FastMath.hh
    FieldDataTypeTraits.hh
    Functors.hh
    GeometricUtilities.hh
    KeyTraits.hh
    NodeListHashMap.hh
<<<<<<< HEAD
=======
    OpenMP_wrapper.hh
>>>>>>> origin/master
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
    intpow2.hh
    iterateIdealH.hh
    lineSegmentIntersections.hh
    log2.hh
    mortonOrderIndices.hh
    mpiUtilities.hh
    newtonRaphson.hh
    nodeBoundingBoxes.hh
    nodeBoundingBoxesInline.hh
    overlayRemapFields.hh
    packElement.hh
    peanoHilbertOrderIndices.hh
    planarReflectingOperator.hh
    pointDistances.hh
    pointInPolygon.hh
    pointInPolyhedron.hh
    pointOnPolygon.hh
    pointOnPolyhedron.hh
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

<<<<<<< HEAD
set(Utilities_python_targets
    fitspline.py
    )

=======
>>>>>>> origin/master
blt_add_library(NAME        Spheral_Utilities
                SOURCES     ${Utilities_sources}
                HEADERS     ${Utilities_headers}
                DEPENDS_ON  ${spheral_depends}
<<<<<<< HEAD
                SHARED      TRUE
=======
>>>>>>> origin/master
                )

install(TARGETS              Spheral_Utilities
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Utilities_headers}
        DESTINATION include/Utilities
        )
<<<<<<< HEAD

install(FILES       ${Utilities_python_targets}
        DESTINATION .
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Utilities)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
=======
>>>>>>> origin/master
