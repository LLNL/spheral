text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/RectangularFinitePlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RectangularFinitePlane< Dim< %(ndim)s > >;
}
"""
