text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidNodeList.cc"

namespace Spheral {
  template class SolidNodeList< Dim< %(ndim)s > >;
}
"""
