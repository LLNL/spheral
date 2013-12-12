//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PointPotential.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class PointPotential< Dim<1> >;
    template class PointPotential< Dim<2> >;
    template class PointPotential< Dim<3> >;
  }
}
