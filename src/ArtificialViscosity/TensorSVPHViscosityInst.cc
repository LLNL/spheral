//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorSVPHViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorSVPHViscosity< Dim<1> >;
    template class TensorSVPHViscosity< Dim<2> >;
    template class TensorSVPHViscosity< Dim<3> >;
  }
}
