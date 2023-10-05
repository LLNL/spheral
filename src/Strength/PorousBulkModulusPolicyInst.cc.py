text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousBulkModulusPolicy.cc"

namespace Spheral {
  template class PorousBulkModulusPolicy<Dim< %(ndim)s > >;
}
"""
