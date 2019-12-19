text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/SumVoronoiMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SumVoronoiMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
