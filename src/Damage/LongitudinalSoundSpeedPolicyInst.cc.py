text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/LongitudinalSoundSpeedPolicy.cc"

namespace Spheral {
  template class LongitudinalSoundSpeedPolicy<Dim< %(ndim)s > >;
}
"""
