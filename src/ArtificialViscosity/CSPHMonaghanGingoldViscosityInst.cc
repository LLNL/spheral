//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CSPHMonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class CSPHMonaghanGingoldViscosity< Dim<1> >;
    template class CSPHMonaghanGingoldViscosity< Dim<2> >;
    template class CSPHMonaghanGingoldViscosity< Dim<3> >;
  }
}
