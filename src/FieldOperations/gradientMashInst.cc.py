text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/gradientMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== gradientMash() ==============================
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>::GradientType> 
gradientMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                     const TableKernel< Dim< %(ndim)s > >& kernel);
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>::GradientType> 
gradientMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                     const TableKernel< Dim< %(ndim)s > >& kernel);

//============================== gradientMash() ==============================
template 
FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>::GradientType> 
gradientMash2<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weightDensity,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                      const TableKernel< Dim< %(ndim)s > >& kernel);

}
"""
