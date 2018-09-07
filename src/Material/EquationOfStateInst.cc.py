text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "EquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class EquationOfState< Dim< %(ndim)s >  >;
}

"""
