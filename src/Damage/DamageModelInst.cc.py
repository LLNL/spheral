text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "DamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DamageModel<Dim< %(ndim)s > >;
}
"""
