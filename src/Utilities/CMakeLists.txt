set(Utility_inst
    Bessel
    boundingBox
    buildMessage
    clipFacetedVolume
    coarsenBinnedValues
    coarsenBinnedValuesInst
    #computeConvexHull
    computeShepardsInterpolation
    DBC
    globalBoundingVolumes
    integrateThroughMeshAlongSegment
    iterateIdealH
    KeyTraits
    lineSegmentIntersections
    medianPosition
    mortonOrderIndices
    nodeBoundingBoxes
    nodeOrdering
    numberDensity
    overlayRemapFields
    peanoHilbertOrderIndices
    PeanoHilbertTransform
    pointInPolygon
    pointInPolyhedron
    pointOnPolygon
    pointOnPolyhedron
    Process
    RedistributionRegistrar
    refinePolyhedron
    RegisterMPIDataTypes
    segmentIntersectEdges
    SpheralTimers
    SpheralTimersAlternate
    Timer
    vectorBoundingBoxInst
    )

set(Utilitiy_sources )
instantiate(Utility_inst Utility_sources)

set(Utility_headers
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

set(Utility_includes
  ${SPHERAL_SOURCE_DIR}
  ${EIGEN_LIB_LOCATION}
  )

blt_add_library(NAME	    Spheral_Utility
		HEADERS     ${Utility_headers}
		SOURCES     ${Utility_sources}
		DEPENDS_ON  ${spheral_depends}
		INCLUDES    ${Utility_includes}
		)

install(TARGETS		    Spheral_Utility
	EXPORT		    ${PROJECT_NAME}-targets
	DESTINATION	    lib
	INCLUDES DESTINATION include
	)

install(FILES	    ${Utility_headers}
	DESTINATION include/Utility
	)	
