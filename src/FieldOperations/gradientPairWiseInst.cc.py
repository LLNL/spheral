text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/gradientPairWise.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== gradient() ==============================
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>::GradientType> 
gradientPairWise<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                         const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>::GradientType> 
gradientPairWise<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                         const TableKernel< Dim< %(ndim)s > >& kernel);

}
"""
