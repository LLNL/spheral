set(Neighbor_sources2d
    ConnectivityMapInst.cc.2.cc
    GridCellIndexInst.cc.2.cc
    GridCellPlaneInst.cc.2.cc
    NeighborInst.cc.2.cc
    NestedGridNeighborInst.cc.2.cc
    TreeNeighborInst.cc.2.cc
    )

set(Neighbor_sources3d
    ConnectivityMapInst.cc.3.cc
    GridCellIndexInst.cc.3.cc
    GridCellPlaneInst.cc.3.cc
    NeighborInst.cc.3.cc
    NestedGridNeighborInst.cc.3.cc
    TreeNeighborInst.cc.3.cc
    )

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
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Neighbor
                SOURCES     ${Neighbor_sources2d} ${Neighbor_sources3d}
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
