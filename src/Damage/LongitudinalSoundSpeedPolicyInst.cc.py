text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "LongitudinalSoundSpeedPolicy.cc"

namespace Spheral {
  template class LongitudinalSoundSpeedPolicy<Dim< %(ndim)s > >;
}
"""
