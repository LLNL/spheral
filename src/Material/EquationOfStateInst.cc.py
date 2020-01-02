text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/EquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class EquationOfState< Dim< %(ndim)s >  >;
}

"""
