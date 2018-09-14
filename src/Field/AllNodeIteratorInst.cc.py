text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/AllNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class AllNodeIterator< Dim< %(ndim)s > >;
}
"""
