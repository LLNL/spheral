text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousSoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousSoundSpeedPolicy<Dim< %(ndim)s > >;
}

"""
