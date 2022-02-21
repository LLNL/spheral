text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/HerzianDEM.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HerzianDEM< Dim< %(ndim)s > >;
}
"""
