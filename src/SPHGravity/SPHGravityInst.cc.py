text = """
//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "SPHGravity.cc"

namespace Spheral {
  template class SPHGravity<Dim< %(ndim)s > >;
} // end namespace Spheral

"""
