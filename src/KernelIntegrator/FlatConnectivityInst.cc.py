text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "FlatConnectivity.cc"

namespace Spheral {
  template class FlatConnectivity<Dim<%(ndim)s>>;
}
"""
