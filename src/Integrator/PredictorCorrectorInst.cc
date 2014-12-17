//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PredictorCorrector.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class PredictorCorrector< Dim<1> >;
    template class PredictorCorrector< Dim<2> >;
    template class PredictorCorrector< Dim<3> >;
  }
}
