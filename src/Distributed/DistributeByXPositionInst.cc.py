text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DistributeByXPosition.cc"

namespace Spheral {
  template class DistributeByXPosition< Dim< %(ndim)s > >;
}

"""
