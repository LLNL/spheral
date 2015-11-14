text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "computeCRKSPHSumMassDensity.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeCRKSPHSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                              const TableKernel<Dim< %(ndim)s > >&, 
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                              const CRKOrder correctionOrder,
                                              FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
  }
}
"""
