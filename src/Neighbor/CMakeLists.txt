blt_add_library( 
  NAME Neighbor_spheral
  HEADERS ConnectivityMap.hh 
	  GridCellIndex.hh
	  GridCellPlane.hh
	  Neighbor.hh
	  NestedGridNeighbor.hh
	  NodePairList.hh
	  TreeNeighbor.hh
  SOURCES ConnectivityMap.cc
	  GridCellIndex.cc
	  GridCellPlane.cc
	  Neighbor.cc
	  NestedGridNeighbor.cc
	  NodePairList.cc
	  TreeNeighbor.cc
  OBJECT FALSE
  DEPENDS_ON Field_spheral)
