text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/RefineNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RefineNodeIterator< Dim< %(ndim)s > >;
}
"""
