//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TaylorSPHHydroBase.cc"

namespace Spheral {
  namespace TaylorSPHSpace {
    template class TaylorSPHHydroBase< Dim<1> >;
    template class TaylorSPHHydroBase< Dim<2> >;
    template class TaylorSPHHydroBase< Dim<3> >;
  }
}
