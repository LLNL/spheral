text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/gradientCRKSPH.cc"
#include "Geometry/Dimension.hh"
#include "SPH/NodeCoupling.hh"

namespace Spheral {

template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>::GradientType> 
gradientCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& gradA,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>& gradC,
                                       const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                       const CRKOrder correctionOrder,
                                       const TableKernel< Dim< %(ndim)s > >& kernel,
                                       const NodeCoupling& nodeCoupling);
  template 
  FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>::GradientType> 
  gradientCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& gradA,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>& gradC,
                                         const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                         const CRKOrder correctionOrder,
                                         const TableKernel< Dim< %(ndim)s > >& kernel,
                                         const NodeCoupling& nodeCoupling);
}

"""
