text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/InternalNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InternalNodeIterator< Dim< %(ndim)s > >;
}
"""
