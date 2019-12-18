text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GenericHydro.cc"

namespace Spheral {
  template class GenericHydro< Dim< %(ndim)s > >;
}
"""
