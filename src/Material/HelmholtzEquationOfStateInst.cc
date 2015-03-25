//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "HelmholtzEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
    namespace Material {
        template class HelmholtzEquationOfState<Dim<1> >;
        template class HelmholtzEquationOfState<Dim<2> >;
        template class HelmholtzEquationOfState<Dim<3> >;
    }
}