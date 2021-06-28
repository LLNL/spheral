text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/computeParticleRadius.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void computeParticleRadius(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                            FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
