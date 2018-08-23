text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PorousStrengthModel.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class PorousStrengthModel<Dim< %(ndim)s > >;
  }
}
"""
