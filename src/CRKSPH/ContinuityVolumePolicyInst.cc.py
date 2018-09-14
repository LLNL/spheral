text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/ContinuityVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ContinuityVolumePolicy<Dim< %(ndim)s > >;
}

"""
