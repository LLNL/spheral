text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/NodeList.cc"

namespace Spheral {
  template class NodeList< Dim< %(ndim)s > >;
}
"""
