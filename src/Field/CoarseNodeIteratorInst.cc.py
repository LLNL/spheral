text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/CoarseNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CoarseNodeIterator< Dim< %(ndim)s > >;
}
"""
