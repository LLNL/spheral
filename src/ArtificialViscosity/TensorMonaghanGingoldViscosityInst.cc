//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorMonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorMonaghanGingoldViscosity< Dim<1> >;
    template class TensorMonaghanGingoldViscosity< Dim<2> >;
    template class TensorMonaghanGingoldViscosity< Dim<3> >;
  }
}
