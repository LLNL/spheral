text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/YoungsModulusPolicy.cc"

namespace Spheral {
  template class YoungsModulusPolicy<Dim< %(ndim)s > >;
}

"""
