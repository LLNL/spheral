text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "interpolateCRKSPH.cc"
#include "Geometry/Dimension.hh"
#include "SPH/NodeCoupling.hh"

namespace Spheral {

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
interpolateCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                          const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                          const CRKOrder correctionOrder,
                                          const TableKernel< Dim< %(ndim)s > >& kernel,
                                          const NodeCoupling& nodeCoupling);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
interpolateCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                          const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                          const CRKOrder correctionOrder,
                                          const TableKernel< Dim< %(ndim)s > >& kernel,
                                          const NodeCoupling& nodeCoupling);

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
interpolateCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                          const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                          const CRKOrder correctionOrder,
                                          const TableKernel< Dim< %(ndim)s > >& kernel,
                                          const NodeCoupling& nodeCoupling);

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
interpolateCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                             const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                             const CRKOrder correctionOrder,
                                             const TableKernel< Dim< %(ndim)s > >& kernel,
                                             const NodeCoupling& nodeCoupling);

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor> 
interpolateCRKSPH<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>& fieldList,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                                   const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                                   const CRKOrder correctionOrder,
                                                   const TableKernel< Dim< %(ndim)s > >& kernel,
                                                   const NodeCoupling& nodeCoupling);

}

"""
