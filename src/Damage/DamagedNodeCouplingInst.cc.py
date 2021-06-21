text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Damage/DamagedNodeCoupling.cc"

namespace Spheral {
  template class DamagedNodeCoupling< Dim< %(ndim)s > >;
}
"""
