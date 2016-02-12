text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MasterNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MasterNodeIterator< Dim< %(ndim)s > >;
}
"""
