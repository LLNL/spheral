text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/TensorStrainPolicy.cc"

namespace Spheral {
  template class Spheral::TensorStrainPolicy<Dim< %(ndim)s > >;
}
"""
