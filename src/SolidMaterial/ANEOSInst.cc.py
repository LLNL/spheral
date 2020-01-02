text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/ANEOS.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ANEOS<Dim< %(ndim)s > >;
}
"""
