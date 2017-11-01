text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList.cc"

namespace Spheral {
  namespace NodeSpace {
    template class NodeList< Dim< %(ndim)s > >;
  }
}
"""
