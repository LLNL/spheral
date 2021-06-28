text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Damage/ThreePointDamagedNodeCoupling.cc"

namespace Spheral {
  template class ThreePointDamagedNodeCoupling< Dim< %(ndim)s > >;
}
"""
