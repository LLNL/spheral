text = """
//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "SPHGravity/SPHGravity.cc"

namespace Spheral {
  template class SPHGravity<Dim< %(ndim)s > >;
} // end namespace Spheral

"""
