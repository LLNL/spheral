text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TensorCRKSPHViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class TensorCRKSPHViscosity< Dim< %(ndim)s > >;
}
"""
