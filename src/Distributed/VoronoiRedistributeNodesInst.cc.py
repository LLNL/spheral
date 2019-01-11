text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "VoronoiRedistributeNodes.cc"

namespace Spheral {
  template class VoronoiRedistributeNodes< Dim< %(ndim)s > >;
}
"""
