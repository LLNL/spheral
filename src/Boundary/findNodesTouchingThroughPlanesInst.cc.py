text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/findNodesTouchingThroughPlanes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template
  std::vector<int> 
  findNodesTouchingThroughPlanes(const NodeList<Dim<%(ndim)s>>& nodeList,
                                 const GeomPlane<Dim<%(ndim)s>>& enterPlane,
                                 const GeomPlane<Dim<%(ndim)s>>& exitPlane,
                                 const double hmultiplier);
                         
}
"""
