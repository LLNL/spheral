text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "splatFieldsMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== splatFieldsMash() ==============================
    template 
    FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
    splatFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const TableKernel< Dim< %(ndim)s > >& kernel,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

    template 
    FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
    splatFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const TableKernel< Dim< %(ndim)s > >& kernel,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

    template 
    FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
    splatFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const TableKernel< Dim< %(ndim)s > >& kernel,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

    template 
    FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
    splatFieldsMash<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                               const TableKernel< Dim< %(ndim)s > >& kernel,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

  }
}
"""
