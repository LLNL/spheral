text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConstantAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantAcceleration< Dim< %(ndim)s > >;
}
"""
