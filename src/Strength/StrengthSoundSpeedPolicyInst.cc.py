text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrengthSoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class StrengthSoundSpeedPolicy<Dim< %(ndim)s > >;
}

"""
