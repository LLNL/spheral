text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/ArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialViscosity<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>;
  template class ArtificialViscosity<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>;
}
"""
