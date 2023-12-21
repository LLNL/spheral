text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorosityModel.cc"

namespace Spheral {
  template class PorosityModel<Dim<%(ndim)s>>;
}
"""
