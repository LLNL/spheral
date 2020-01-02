text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/GammaPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GammaPolicy<Dim< %(ndim)s > >;
}

"""
