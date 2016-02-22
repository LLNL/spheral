text = """
//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "SPHGravity.cc"

namespace Spheral {
  namespace GravitySpace {

    template class SPHGravity<Dim< %(ndim)s > >;

  } // end namespace GravitySpace
} // end namespace Spheral

"""
