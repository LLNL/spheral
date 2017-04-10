text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GhostNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GhostNodeIterator< Dim< %(ndim)s > >;
}
"""
