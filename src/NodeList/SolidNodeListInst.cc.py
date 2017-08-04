text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidNodeList.cc"

namespace Spheral {
  namespace NodeSpace {
    template class SolidNodeList< Dim< %(ndim)s > >;
  }
}
"""
