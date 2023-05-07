text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/InfinitePlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InfinitePlane< Dim< %(ndim)s > >;
}
"""
