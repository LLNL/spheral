//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SPHSpace {
    template class SPHHydroBase< Dim<1> >;
    template class SPHHydroBase< Dim<2> >;
    template class SPHHydroBase< Dim<3> >;
  }
}
