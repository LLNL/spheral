text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GridCellPlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NeighborSpace {
    template class GridCellPlane< Dim< %(ndim)s > >;
  }
}
"""
