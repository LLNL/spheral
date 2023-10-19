text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousSoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousSoundSpeedPolicy<Dim< %(ndim)s > >;
}

"""
