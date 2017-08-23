text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computeFragmentField.cc"

namespace Spheral {
  template Field<Dim< %(ndim)s >, int> computeFragmentField(const NodeList<Dim< %(ndim)s > >&, const double, const Field<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&, const double, const bool);
}
"""
