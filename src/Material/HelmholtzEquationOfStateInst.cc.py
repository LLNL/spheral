text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "HelmholtzEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
    namespace Material {
        template class HelmholtzEquationOfState< Dim< %(ndim)s > >;
    }
}
"""
