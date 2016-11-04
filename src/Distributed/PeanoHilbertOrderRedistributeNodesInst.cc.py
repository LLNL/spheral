text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PeanoHilbertOrderRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class PeanoHilbertOrderRedistributeNodes< Dim< %(ndim)s > >;
  }
}

"""
