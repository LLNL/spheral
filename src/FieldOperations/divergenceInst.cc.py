text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/divergence.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== divergence() ==============================
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>::DivergenceType> 
divergence<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                   const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>::DivergenceType> 
divergence<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                   const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>::DivergenceType> 
divergence<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                      const TableKernel< Dim< %(ndim)s > >& kernel);

}
"""
