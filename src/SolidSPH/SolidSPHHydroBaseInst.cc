//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidSPHSpace {
    template class SolidSPHHydroBase< Dim<1> >;
    template class SolidSPHHydroBase< Dim<2> >;
    template class SolidSPHHydroBase< Dim<3> >;
  }
}
