text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "JohnsonCookFailureStrainPolicy.cc"

namespace Spheral {
  template class Spheral::JohnsonCookFailureStrainPolicy<Dim< %(ndim)s > >;
}
"""
