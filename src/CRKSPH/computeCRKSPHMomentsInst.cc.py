text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/computeCRKSPHMoments.cc"

namespace Spheral {

template void computeCRKSPHMoments(const ConnectivityMap< Dim< %(ndim)s > >& connectivityMap,
                                   const TableKernel< Dim< %(ndim)s > >& W,
                                   const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Scalar>& weight,
                                   const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& position,
                                   const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::SymTensor>& H,
                                   const RKOrder correctionOrder,
                                   const NodeCoupling& nodeCoupling,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Scalar>& m0,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& m1,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::SymTensor>& m2,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::ThirdRankTensor>& m3,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::FourthRankTensor>& m4,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& gradm0,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Tensor>& gradm1,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::ThirdRankTensor>& gradm2,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::FourthRankTensor>& gradm3,
                                   FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::FifthRankTensor>& gradm4);

}

"""
