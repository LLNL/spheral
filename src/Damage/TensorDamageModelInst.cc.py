text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "TensorDamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class TensorDamageModel<Dim< %(ndim)s > >;
}
"""
