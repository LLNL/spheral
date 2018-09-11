text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Hydro/SecondMomentHourglassControl.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class SecondMomentHourglassControl< Dim< %(ndim)s > >;
  }
}
"""
