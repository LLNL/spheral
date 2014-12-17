//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorCSPHViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorCSPHViscosity< Dim<1> >;
    template class TensorCSPHViscosity< Dim<2> >;
    template class TensorCSPHViscosity< Dim<3> >;
  }
}
