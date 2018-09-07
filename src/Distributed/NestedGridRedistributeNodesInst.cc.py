text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NestedGridRedistributeNodes.cc"

namespace Spheral {
  template class NestedGridRedistributeNodes< Dim< %(ndim)s > >;
}

"""
