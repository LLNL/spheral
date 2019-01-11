text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/MasterNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MasterNodeIterator< Dim< %(ndim)s > >;
}
"""
