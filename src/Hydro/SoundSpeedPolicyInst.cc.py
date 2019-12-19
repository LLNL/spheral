text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/SoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SoundSpeedPolicy<Dim< %(ndim)s > >;
}

"""
