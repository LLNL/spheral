text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/Policies/ALEPositionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ALEPositionPolicy<Dim< %(ndim)s > >;
}
"""
