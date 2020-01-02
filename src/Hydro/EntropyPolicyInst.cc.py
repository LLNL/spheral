text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/EntropyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class EntropyPolicy<Dim< %(ndim)s > >;
}

"""
