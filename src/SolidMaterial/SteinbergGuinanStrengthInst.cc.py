text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/SteinbergGuinanStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SteinbergGuinanStrength<Dim< %(ndim)s > >;
}
"""
