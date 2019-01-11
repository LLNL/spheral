text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList.cc"

namespace Spheral {
  template class NodeList< Dim< %(ndim)s > >;
}
"""
