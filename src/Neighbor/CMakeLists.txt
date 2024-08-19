include_directories(.)

set(Neighbor_sources 
    NodePairList.cc)

# Add explicit instantiations
list(APPEND Neighbor_sources
    ConnectivityMapInst.cc
    GridCellIndexInst.cc
    GridCellPlaneInst.cc
    NeighborInst.cc
    NestedGridNeighborInst.cc
    TreeNeighborInst.cc
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
    NodePairList.hh
    TreeNeighbor.hh
    )

spheral_add_obj_library(Neighbor SPHERAL_OBJ_LIBS)

