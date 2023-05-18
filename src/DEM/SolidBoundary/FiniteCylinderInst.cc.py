text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/FiniteCylinder.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FiniteCylinder< Dim< %(ndim)s > >;
}
"""
