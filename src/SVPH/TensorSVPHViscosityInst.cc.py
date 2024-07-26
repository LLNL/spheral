text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SVPH/TensorSVPHViscosity.cc"

namespace Spheral {
  template class TensorSVPHViscosity< Dim< %(ndim)s > >;
}
"""
