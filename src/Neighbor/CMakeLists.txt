set(Neighbor_inst
    ConnectivityMap
    GridCellIndex
    GridCellPlane
    Neighbor
    NestedGridNeighbor
    NodePairList
    TreeNeighbor
    )

set(Neighbor_sources )
instantiate(Neighbor_inst Neighbor_sources)

set(Neighbor_headers
    ConnectivityMap.hh
    ConnectivityMapInline.hh
    GridCellIndex.hh
    GridCellIndexBase.hh
    GridCellIndexInline.hh
    GridCellPlane.hh
    GridCellPlaneInline.hh
    Neighbor.hh
    NeighborInline.hh
    NestedGridNeighbor.hh
    NestedGridNeighborInline.hh
    NodePairList.hh
    TreeNeighbor.hh
    )

set(Neighbor_includes
  ${SPHERAL_SOURCE_DIR}
  ${EIGEN_LIB_LOCATION}
  )
 
blt_add_library(NAME	    Spheral_Neighbor
		HEADERS     ${Neighbor_headers}
		SOURCES     ${Neighbor_sources}
		DEPENDS_ON  ${spheral_depends}
		INCLUDES    ${Neighbor_includes}
		)

install(TARGETS		    Spheral_Neighbor
	EXPORT		    ${PROJECT_NAME}-targets
	DESTINATION	    lib
	INCLUDES DESTINATION include
	)

install(FILES	    ${Neighbor_headers}
	DESTINATION include/Neighbor
	)
