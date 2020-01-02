text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/VolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VolumePolicy<Dim< %(ndim)s > >;
}

"""
