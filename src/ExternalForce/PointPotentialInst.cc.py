text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PointPotential.cc"

namespace Spheral {
  template class PointPotential< Dim< %(ndim)s > >;
}
"""
