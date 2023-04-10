text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/SolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidBoundary< Dim< %(ndim)s > >;
}
"""
