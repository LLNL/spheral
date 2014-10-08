//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SVPHFacetedHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SVPHSpace {
    template class SVPHFacetedHydroBase< Dim<1> >;
    template class SVPHFacetedHydroBase< Dim<2> >;
    template class SVPHFacetedHydroBase< Dim<3> >;
  }
}
