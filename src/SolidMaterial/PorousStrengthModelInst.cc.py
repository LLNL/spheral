text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PorousStrengthModel.cc"

namespace Spheral {
  template class PorousStrengthModel<Dim< %(ndim)s > >;
}
"""
