text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NestedGridRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class NestedGridRedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
