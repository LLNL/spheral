text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/StrengthModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class StrengthModel<Dim< %(ndim)s > >;
}
"""
