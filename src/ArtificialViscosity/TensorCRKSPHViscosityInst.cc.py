text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TensorCRKSPHViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorCRKSPHViscosity< Dim< %(ndim)s > >;
  }
}
"""
