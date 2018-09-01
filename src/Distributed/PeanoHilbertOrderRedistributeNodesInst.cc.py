text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PeanoHilbertOrderRedistributeNodes.cc"

namespace Spheral {
  template class PeanoHilbertOrderRedistributeNodes< Dim< %(ndim)s > >;
}

"""
