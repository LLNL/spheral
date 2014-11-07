//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class Integrator< Dim<1> >;
    template class Integrator< Dim<2> >;
    template class Integrator< Dim<3> >;
  }
}
