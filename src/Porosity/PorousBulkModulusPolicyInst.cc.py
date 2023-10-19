text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousBulkModulusPolicy.cc"

namespace Spheral {
  template class PorousBulkModulusPolicy<Dim< %(ndim)s > >;
}
"""
