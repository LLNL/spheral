text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/SolidNodeList.cc"

namespace Spheral {
  template class SolidNodeList< Dim< %(ndim)s > >;
}
"""
