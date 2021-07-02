text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Damage/DamageGradientNodeCoupling.cc"

namespace Spheral {
  template class DamageGradientNodeCoupling< Dim< %(ndim)s > >;
}
"""
