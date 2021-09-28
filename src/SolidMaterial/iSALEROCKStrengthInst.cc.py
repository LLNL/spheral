text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/iSALEROCKStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class iSALEROCKStrength<Dim< %(ndim)s > >;
}
"""
