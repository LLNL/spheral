text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPHMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
