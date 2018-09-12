text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/divergenceMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== divergenceMash() ==============================
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>::DivergenceType> 
divergenceMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                       const TableKernel< Dim< %(ndim)s > >& kernel);

template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>::DivergenceType> 
divergenceMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                       const TableKernel< Dim< %(ndim)s > >& kernel);

}
"""
