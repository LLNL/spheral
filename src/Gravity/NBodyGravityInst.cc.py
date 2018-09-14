text = """
//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "Gravity/NBodyGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace GravitySpace {

    template class NBodyGravity<Dim< %(ndim)s > >;

  } // end namespace GravitySpace
} // end namespace Spheral

"""
