text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "BilinearIndex.cc"

namespace Spheral {
  template class BilinearIndex<Dim<%(ndim)s>>;
}
"""
