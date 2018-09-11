text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Hydro/ThirdMomentHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class ThirdMomentHourglassControl< Dim< %(ndim)s > >;
  }
}
"""
