text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorSVPHViscosity.cc"

namespace Spheral {
  template class TensorSVPHViscosity< Dim< %(ndim)s > >;
}
"""
