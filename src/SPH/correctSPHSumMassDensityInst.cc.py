text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/correctSPHSumMassDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void correctSPHSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                         const TableKernel<Dim< %(ndim)s > >&, 
                                         const bool,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                         FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
