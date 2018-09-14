text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/MortonOrderRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class MortonOrderRedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
