text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/IvanoviSALEDamagePolicy.cc"

namespace Spheral {
  template class IvanoviSALEDamagePolicy<Dim< %(ndim)s > >;
}
"""
