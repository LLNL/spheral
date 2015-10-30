text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "InternalNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InternalNodeIterator< Dim< %(ndim)s > >;
}
"""
