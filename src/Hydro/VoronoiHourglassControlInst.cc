//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "VoronoiHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class VoronoiHourglassControl< Dim<1> >;
    template class VoronoiHourglassControl< Dim<2> >;
    template class VoronoiHourglassControl< Dim<3> >;
  }
}
