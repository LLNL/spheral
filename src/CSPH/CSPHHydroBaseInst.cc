//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CSPHHydroBase.cc"

namespace Spheral {
  namespace CSPHSpace {
    template class CSPHHydroBase< Dim<1> >;
    template class CSPHHydroBase< Dim<2> >;
    template class CSPHHydroBase< Dim<3> >;
  }
}
