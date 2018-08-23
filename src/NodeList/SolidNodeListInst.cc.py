text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/SolidNodeList.cc"

namespace Spheral {
  namespace NodeSpace {
    template class SolidNodeList< Dim< %(ndim)s > >;
  }
}
"""
