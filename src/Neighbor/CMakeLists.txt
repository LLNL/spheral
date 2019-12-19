set(Neighbor_inst
    ConnectivityMap
    GridCellIndex
    GridCellPlane
    Neighbor
    NestedGridNeighbor
    TreeNeighbor
   )


set(Neighbor_sources
    NodePairList.cc
   )

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

blt_add_library(NAME        Spheral_Neighbor
                SOURCES     ${Neighbor_sources}
                HEADERS     ${Neighbor_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Neighbor
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Neighbor_headers}
        DESTINATION include/Neighbor
        )
