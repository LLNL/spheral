text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "VoronoiRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class VoronoiRedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
