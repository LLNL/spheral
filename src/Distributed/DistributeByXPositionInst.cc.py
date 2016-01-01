text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DistributeByXPosition.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class DistributeByXPosition< Dim< %(ndim)s > >;
  }
}

"""
