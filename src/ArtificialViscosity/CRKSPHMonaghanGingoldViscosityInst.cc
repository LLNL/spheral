//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPHMonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class CRKSPHMonaghanGingoldViscosity< Dim<1> >;
    template class CRKSPHMonaghanGingoldViscosity< Dim<2> >;
    template class CRKSPHMonaghanGingoldViscosity< Dim<3> >;
  }
}
