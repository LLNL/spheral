//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "MorrisMonaghanReducingViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class MorrisMonaghanReducingViscosity< Dim<1> >;
    template class MorrisMonaghanReducingViscosity< Dim<2> >;
    template class MorrisMonaghanReducingViscosity< Dim<3> >;
  }
}
