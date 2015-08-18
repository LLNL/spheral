text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPHCorrectionsPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHCorrectionsPolicy<Dim< %(ndim)s > >;
}

"""
