//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScaleBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NodeSpace {
    template class SmoothingScaleBase< Dim<1> >;
    template class SmoothingScaleBase< Dim<2> >;
    template class SmoothingScaleBase< Dim<3> >;
  }
}
