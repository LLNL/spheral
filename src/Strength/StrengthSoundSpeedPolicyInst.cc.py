text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/StrengthSoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class StrengthSoundSpeedPolicy<Dim< %(ndim)s > >;
}

"""
