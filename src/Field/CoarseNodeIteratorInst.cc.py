text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CoarseNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CoarseNodeIterator< Dim< %(ndim)s > >;
}
"""
