text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SumVoronoiMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SumVoronoiMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
