//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ThirdMomentHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class ThirdMomentHourglassControl< Dim<1> >;
    template class ThirdMomentHourglassControl< Dim<2> >;
    template class ThirdMomentHourglassControl< Dim<3> >;
  }
}
