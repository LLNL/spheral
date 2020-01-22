text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/sampleMultipleFieldsMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== sampleFieldsMash() ==============================
template 
FieldListSet< Dim< %(ndim)s > >
sampleMultipleFieldsMash< Dim< %(ndim)s > >(const FieldListSet< Dim< %(ndim)s > >& fieldListSet,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                   const TableKernel< Dim< %(ndim)s > >& kernel,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield);

}
"""
