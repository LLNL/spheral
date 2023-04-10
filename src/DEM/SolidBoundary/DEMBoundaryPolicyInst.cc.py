text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DEM/SolidBoundary/DEMBoundaryPolicy.cc"

namespace Spheral {
  template class DEMBoundaryPolicy<Dim< %(ndim)s >>;
}
"""
