text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Damage/PairMaxDamageNodeCoupling.cc"

namespace Spheral {
  template class PairMaxDamageNodeCoupling< Dim< %(ndim)s > >;
}
"""
