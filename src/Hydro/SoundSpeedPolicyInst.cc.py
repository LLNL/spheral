text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SoundSpeedPolicy<Dim< %(ndim)s > >;
}

"""
