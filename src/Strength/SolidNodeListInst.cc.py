text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidNodeList.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class SolidNodeList< Dim< %(ndim)s > >;
  }
}
"""
