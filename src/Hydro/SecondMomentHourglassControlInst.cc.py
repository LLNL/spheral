text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SecondMomentHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class SecondMomentHourglassControl< Dim< %(ndim)s > >;
  }
}
"""
