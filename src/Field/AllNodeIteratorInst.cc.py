text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "AllNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class AllNodeIterator< Dim< %(ndim)s > >;
}
"""
