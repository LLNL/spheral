text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SVPHMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
