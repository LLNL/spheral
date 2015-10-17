text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHCorrections.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           const int correctionOrder,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);

    template void computeZerothCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&);

    template void computeLinearCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);

    template void computeQuadraticCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);

    template void computeQuadraticCRKSPHCorrectionsMike(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);

    template void computeCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           const NodeCoupling& nodeCoupling,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);
  }
}

"""
