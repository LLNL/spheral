text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PointPotential.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class PointPotential< Dim< %(ndim)s > >;
  }
}
"""
