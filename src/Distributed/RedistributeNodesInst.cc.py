text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class RedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
