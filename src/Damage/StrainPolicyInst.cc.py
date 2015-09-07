text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "StrainPolicy.cc"

namespace Spheral {
  template class Spheral::StrainPolicy<Dim< %(ndim)s > >;
}
"""
