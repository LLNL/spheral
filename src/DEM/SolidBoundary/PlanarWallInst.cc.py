text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/PlanarWall.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PlanarWall< Dim< %(ndim)s > >;
}
"""
