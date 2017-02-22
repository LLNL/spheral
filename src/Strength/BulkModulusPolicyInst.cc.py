text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "BulkModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class BulkModulusPolicy<Dim< %(ndim)s > >;
}
"""
