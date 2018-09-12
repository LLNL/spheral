text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ExternalForce/ConstantAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantAcceleration< Dim< %(ndim)s > >;
}
"""
