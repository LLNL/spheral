text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/smoothFields.cc"
#include "Geometry/Dimension.hh"

//============================== smoothFields() ==============================
namespace Spheral {

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
smoothFields<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                                       const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
smoothFields<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                                       const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
smoothFields<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                                       const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
smoothFields<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                                          const TableKernel< Dim< %(ndim)s > >& kernel);

}
"""
