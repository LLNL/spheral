text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/TensorDamagePolicy.cc"

namespace Spheral {
  template class TensorDamagePolicy<Dim< %(ndim)s > >;
}
"""
