text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SVPHCorrectionsPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHCorrectionsPolicy<Dim< %(ndim)s > >;
}

"""
