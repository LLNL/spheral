text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/HVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HVolumePolicy<Dim< %(ndim)s > >;
}

"""
