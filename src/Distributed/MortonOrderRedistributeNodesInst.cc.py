text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "MortonOrderRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class MortonOrderRedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
