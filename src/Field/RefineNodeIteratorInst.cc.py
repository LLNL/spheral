text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RefineNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RefineNodeIterator< Dim< %(ndim)s > >;
}
"""
