//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SVPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SVPHSpace {
    template class SVPHHydroBase< Dim<1> >;
    template class SVPHHydroBase< Dim<2> >;
    template class SVPHHydroBase< Dim<3> >;
  }
}
