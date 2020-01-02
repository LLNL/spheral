text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/ContinuityVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ContinuityVolumePolicy<Dim< %(ndim)s > >;
}

"""
