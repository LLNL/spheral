text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PorousStrengthModel.cc"

namespace Spheral {
  template class PorousStrengthModel<Dim< %(ndim)s > >;
}
"""
