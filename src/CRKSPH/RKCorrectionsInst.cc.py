text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RKCorrections.cc"

namespace Spheral {
  template class RKCorrections<Dim<%(ndim)s>>;
}
"""
