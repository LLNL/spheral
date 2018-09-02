text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SecondMomentHourglassControl.cc"

namespace Spheral {
  template class SecondMomentHourglassControl< Dim< %(ndim)s > >;
}
"""
