text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHCorrections.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeCRKSPHCorrections(const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& m1,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& m2,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& m3,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& m4,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Tensor>& gradm1,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& gradm2,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& gradm3,
                                           const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FifthRankTensor>& gradm4,
                                           const CRKOrder correctionOrder,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);

    template void computeZerothCRKSPHCorrections(const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                                 const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&);

    template void computeLinearCRKSPHCorrections(const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                                 const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& m1,
                                                 const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& m2,
                                                 const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                                 const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Tensor>& gradm1,
                                                 const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& gradm2,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);

    template void computeQuadraticCRKSPHCorrections(const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& m1,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& m2,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& m3,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& m4,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Tensor>& gradm1,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& gradm2,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& gradm3,
                                                    const FieldSpace::FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FifthRankTensor>& gradm4,
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
