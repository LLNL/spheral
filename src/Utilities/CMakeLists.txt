set(Utilities_inst
    #Bessel
    boundingBox
    #buildMessage
    #coarsenBinnedValues
    #computeConvexHull
    computeShepardsInterpolation
    globalBoundingVolumes
    integrateThroughMeshAlongSegment
    iterateIdealH
    medianPosition
    mortonOrderIndices
    nodeBoundingBoxes
    nodeOrdering
    numberDensity
    #overlayRemapFields
    peanoHilbertOrderIndices
    #SpheralTimers
    #SpheralTimersAlternate
    #Timer
    )

set(Utilitiy_sources
    DBC
    Process
    KeyTraits
    PeanoHilbertTransform
    RegisterMPIDataTypes
    RedistributionRegistrar
    vectorBoundingBoxInst
    lineSegmentIntersections
    segmentIntersectEdges
    pointInPolygon
    pointInPolyhedron
    pointOnPolygon
    pointOnPolyhedron
    refinePolyhedron
    coarsenBinnedValuesInst
    clipFacetedVolume
  )

instantiate(Utilities_inst Utilities_sources)

set(Utilities_headers
    allReduce.hh
    Bessel.hh
    bisectRoot.hh
    bisectSearch.hh
    boundingBox.hh
    boundPointWithinBox.hh
    clipFacetedVolume.hh
    coarsenBinnedValues.hh
    comparisons.hh
    #computeConvexHull.hh
    computeShepardsInterpolation.hh
    CounterClockwiseComparator.hh
    DataTypeTraits.hh
    DBC.hh
    erff.hh
    fakestl.hh
    FastMath.hh
    FieldDataTypeTraits.hh
    Functors.hh
    GeometricUtilities.hh
    globalBoundingVolumes.hh
    globalNodeIDs.hh
    globalNodeIDsInline.hh
    hashes.hh
    initializeTau.hh
    integrateThroughMeshAlongSegment.hh
    intpow2.hh
    iterateIdealH.hh
    KeyTraits.hh
    lineSegmentIntersections.hh
    log2.hh
    medianPosition.hh
    mortonOrderIndices.hh
    mpiUtilities.hh
    newtonRaphson.hh
    nodeBoundingBoxes.hh
    nodeBoundingBoxesInline.hh
    NodeListHashMap.hh
    nodeOrdering.hh
    numberDensity.hh
    overlayRemapFields.hh
    packElement.hh
    PairComparisons.hh
    peanoHilbertOrderIndices.hh
    PeanoHilbertTransform.hh
    planarReflectingOperator.hh
    pointDistances.hh
    pointInPolygon.hh
    pointInPolyhedron.hh
    pointOnPolygon.hh
    pointOnPolyhedron.hh
    Process.hh
    RedistributionNotification.hh
    RedistributionNotificationHandle.hh
    RedistributionRegistrar.hh
    RedistributionRegistrarInline.hh
    refinePolyhedron.hh
    RegisterMPIDataTypes.hh
    registerWithRedistribution.hh
    removeElements.hh
    rotationMatrix.hh
    SafeIndexMap.hh
    safeInv.hh
    scan.hh
    segmentIntersectEdges.hh
    simpsonsIntegration.hh
    sort_permutation.hh
    SpheralFunctions.hh
    testBoxIntersection.hh
    Timer.hh
    timingUtilities.hh
    Tree.hh
    TreeInline.hh
    )

set(Utilities_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Utilities
  ${EIGEN_LIB_LOCATION}
  )

blt_add_library(NAME        Spheral_Utilities
                HEADERS     ${Utilities_headers}
                SOURCES     ${Utilities_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Utilities_includes}
                )

install(TARGETS             Spheral_Utilities
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Utilities_headers}
        DESTINATION include/Utilities
        )       
