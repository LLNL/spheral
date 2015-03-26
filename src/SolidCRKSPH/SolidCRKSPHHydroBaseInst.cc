//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidCRKSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidCRKSPHSpace {
    template class SolidCRKSPHHydroBase< Dim<1> >;
    template class SolidCRKSPHHydroBase< Dim<2> >;
    template class SolidCRKSPHHydroBase< Dim<3> >;
  }
}
