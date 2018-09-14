text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/JohnsonCookFailureStrainPolicy.cc"

namespace Spheral {
  template class Spheral::JohnsonCookFailureStrainPolicy<Dim< %(ndim)s > >;
}
"""
