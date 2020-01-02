text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/PeanoHilbertOrderRedistributeNodes.cc"

namespace Spheral {
  template class PeanoHilbertOrderRedistributeNodes< Dim< %(ndim)s > >;
}

"""
