//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "VonNeumanViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class VonNeumanViscosity< Dim<1> >;
    template class VonNeumanViscosity< Dim<2> >;
    template class VonNeumanViscosity< Dim<3> >;
  }
}
