text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Hydro/SecondMomentHourglassControl.cc"

namespace Spheral {
  template class SecondMomentHourglassControl< Dim< %(ndim)s > >;
}
"""
