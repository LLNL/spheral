text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorSVPHViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorSVPHViscosity< Dim< %(ndim)s > >;
  }
}
"""
