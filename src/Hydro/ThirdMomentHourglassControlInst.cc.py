text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Hydro/ThirdMomentHourglassControl.cc"

namespace Spheral {
  template class ThirdMomentHourglassControl< Dim< %(ndim)s > >;
}
"""
