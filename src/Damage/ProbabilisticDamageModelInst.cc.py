text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Damage/ProbabilisticDamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ProbabilisticDamageModel<Dim< %(ndim)s > >;
}
"""
