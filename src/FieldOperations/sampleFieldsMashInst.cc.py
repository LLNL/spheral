text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/sampleFieldsMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== sampleFieldsMash() ==============================
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
sampleFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                         const TableKernel< Dim< %(ndim)s > >& kernel,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
sampleFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                         const TableKernel< Dim< %(ndim)s > >& kernel,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
sampleFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                         const TableKernel< Dim< %(ndim)s > >& kernel,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
sampleFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const TableKernel< Dim< %(ndim)s > >& kernel,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

}
"""
