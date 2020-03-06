include_directories(.)
set(Neighbor_inst
    ConnectivityMap
    GridCellIndex
    GridCellPlane
    Neighbor
    NestedGridNeighbor
    TreeNeighbor
   )


set(Neighbor_sources 
    NodePairList.cc)

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
    TreeNeighbor.hh
    NodePairList.hh
    )

spheral_add_cxx_library(Neighbor)

install(TARGETS              Spheral_Neighbor
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Neighbor_headers}
        DESTINATION include/Neighbor
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Neighbor)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
