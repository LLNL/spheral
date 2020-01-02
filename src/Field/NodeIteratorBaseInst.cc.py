text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/NodeIteratorBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NodeIteratorBase< Dim< %(ndim)s > >;
}
"""
