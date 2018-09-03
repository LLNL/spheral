text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ANEOS.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ANEOS<Dim< %(ndim)s > >;
}
"""
