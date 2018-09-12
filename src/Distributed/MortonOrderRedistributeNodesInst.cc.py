text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/MortonOrderRedistributeNodes.cc"

namespace Spheral {
  template class MortonOrderRedistributeNodes< Dim< %(ndim)s > >;
}

"""
