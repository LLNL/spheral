text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ContinuityVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ContinuityVolumePolicy<Dim< %(ndim)s > >;
}

"""
