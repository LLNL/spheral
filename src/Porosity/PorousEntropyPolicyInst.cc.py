text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousEntropyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousEntropyPolicy<Dim< %(ndim)s > >;
}

"""
