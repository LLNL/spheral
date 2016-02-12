text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ThirdMomentHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class ThirdMomentHourglassControl< Dim< %(ndim)s > >;
  }
}
"""
