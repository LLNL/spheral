text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PiecewiseLinearPorousStrengthModel.cc"

namespace Spheral {
  template class PiecewiseLinearPorousStrengthModel<Dim< %(ndim)s > >;
}
"""
