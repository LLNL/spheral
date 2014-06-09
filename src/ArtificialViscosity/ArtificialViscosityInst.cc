//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class ArtificialViscosity< Dim<1> >;
    template class ArtificialViscosity< Dim<2> >;
    template class ArtificialViscosity< Dim<3> >;
  }
}
