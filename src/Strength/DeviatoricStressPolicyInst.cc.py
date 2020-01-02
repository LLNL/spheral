text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/DeviatoricStressPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DeviatoricStressPolicy<Dim< %(ndim)s > >;
}
"""
