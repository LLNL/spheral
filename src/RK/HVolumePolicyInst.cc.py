text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/HVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HVolumePolicy<Dim< %(ndim)s > >;
}

"""
