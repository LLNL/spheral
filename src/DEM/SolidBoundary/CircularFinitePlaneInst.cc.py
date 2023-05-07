text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/CircularFinitePlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CircularFinitePlane< Dim< %(ndim)s > >;
}
"""
