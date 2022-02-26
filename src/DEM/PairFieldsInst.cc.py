text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/PairFields.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PairFields< Dim< %(ndim)s > >;
}
"""
