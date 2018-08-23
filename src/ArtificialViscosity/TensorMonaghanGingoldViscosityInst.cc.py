text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class TensorMonaghanGingoldViscosity< Dim< %(ndim)s > >;
  }
}
"""
