text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorSVPHViscosity.cc"

namespace Spheral {
  template class TensorSVPHViscosity< Dim< %(ndim)s > >;
}
"""
