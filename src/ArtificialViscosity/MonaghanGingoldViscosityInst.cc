//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "MonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class MonaghanGingoldViscosity< Dim<1> >;
    template class MonaghanGingoldViscosity< Dim<2> >;
    template class MonaghanGingoldViscosity< Dim<3> >;
  }
}
