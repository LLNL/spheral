text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/NestedGridRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class NestedGridRedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
