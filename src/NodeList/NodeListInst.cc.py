text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/NodeList.cc"

namespace Spheral {
  namespace NodeSpace {
    template class NodeList< Dim< %(ndim)s > >;
  }
}
"""
