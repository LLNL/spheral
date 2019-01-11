text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GridCellPlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GridCellPlane< Dim< %(ndim)s > >;
}
"""
