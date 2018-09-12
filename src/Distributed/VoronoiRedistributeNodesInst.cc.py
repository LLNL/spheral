text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/VoronoiRedistributeNodes.cc"

namespace Spheral {
  template class VoronoiRedistributeNodes< Dim< %(ndim)s > >;
}
"""
