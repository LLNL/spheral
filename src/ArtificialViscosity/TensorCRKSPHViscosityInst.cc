//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorCRKSPHViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorCRKSPHViscosity< Dim<1> >;
    template class TensorCRKSPHViscosity< Dim<2> >;
    template class TensorCRKSPHViscosity< Dim<3> >;
  }
}
