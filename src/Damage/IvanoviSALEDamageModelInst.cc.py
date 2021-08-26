text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Damage/IvanoviSALEDamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class IvanoviSALEDamageModel<Dim< %(ndim)s > >;
}
"""
