//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SecondMomentHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class SecondMomentHourglassControl< Dim<1> >;
    template class SecondMomentHourglassControl< Dim<2> >;
    // template class SecondMomentHourglassControl< Dim<3> >;
  }
}
