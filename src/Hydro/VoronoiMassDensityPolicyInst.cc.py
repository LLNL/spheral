text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VoronoiMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VoronoiMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
