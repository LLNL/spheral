text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/TensorCRKSPHViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class TensorCRKSPHViscosity< Dim< %(ndim)s > >;
}
"""
