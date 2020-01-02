text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/smoothFieldsMash2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== smoothFieldsMash2() ==============================
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
smoothFieldsMash2<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weightDensity,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                          const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
smoothFieldsMash2<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weightDensity,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                          const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
smoothFieldsMash2<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weightDensity,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                          const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
smoothFieldsMash2<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weightDensity,
                                             const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                             const TableKernel< Dim< %(ndim)s > >& kernel);

}
"""
