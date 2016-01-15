text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorStrainPolicy.cc"

namespace Spheral {
  template class Spheral::TensorStrainPolicy<Dim< %(ndim)s > >;
}
"""
