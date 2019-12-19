text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/BulkModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class BulkModulusPolicy<Dim< %(ndim)s > >;
}
"""
