text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/DampedLinearSpring.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DampedLinearSpring< Dim< %(ndim)s > >;
}
"""
