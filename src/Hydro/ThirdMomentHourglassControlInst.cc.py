text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ThirdMomentHourglassControl.cc"

namespace Spheral {
  template class ThirdMomentHourglassControl< Dim< %(ndim)s > >;
}
"""
