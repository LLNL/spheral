text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PorousStrengthModel.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class PorousStrengthModel<Dim< %(ndim)s > >;
  }
}
"""
