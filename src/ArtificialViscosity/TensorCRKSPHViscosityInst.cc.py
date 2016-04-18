text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorCRKSPHViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorCRKSPHViscosity< Dim< %(ndim)s > >;
  }
}
"""
