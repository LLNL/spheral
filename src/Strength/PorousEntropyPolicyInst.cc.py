text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousEntropyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousEntropyPolicy<Dim< %(ndim)s > >;
}

"""
