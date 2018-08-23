text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/DistributeByXPosition.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class DistributeByXPosition< Dim< %(ndim)s > >;
  }
}

"""
