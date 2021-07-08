text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/ProbabilisticDamagePolicy.cc"

namespace Spheral {
  template class ProbabilisticDamagePolicy<Dim< %(ndim)s > >;
}
"""
