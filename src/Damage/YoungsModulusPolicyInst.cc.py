text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "YoungsModulusPolicy.cc"

namespace Spheral {
  template class YoungsModulusPolicy<Dim< %(ndim)s > >;
}

"""
