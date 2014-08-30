//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "FiniteVolumeViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class FiniteVolumeViscosity< Dim<1> >;
    template class FiniteVolumeViscosity< Dim<2> >;
    template class FiniteVolumeViscosity< Dim<3> >;
  }
}
