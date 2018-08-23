text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/computeHVolumes.cc"

namespace Spheral {
  template void computeHVolumes(const Dim< %(ndim)s >::Scalar nPerh,
                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}

"""
