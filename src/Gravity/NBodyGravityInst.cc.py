text = """
//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "NBodyGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace GravitySpace {

    template class NBodyGravity<Dim< %(ndim)s > >;

  } // end namespace GravitySpace
} // end namespace Spheral

"""
