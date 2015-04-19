//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPHHydroBase.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template class CRKSPHHydroBase< Dim<1> >;
    template class CRKSPHHydroBase< Dim<2> >;
    template class CRKSPHHydroBase< Dim<3> >;
  }
}
