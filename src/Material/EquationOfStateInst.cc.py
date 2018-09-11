text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/EquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class EquationOfState< Dim< %(ndim)s >  >;
  }
}

"""
