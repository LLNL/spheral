text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/sampleMultipleFields2LatticeMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//============================== sampleFieldsMash() ==============================
template 
void
sampleMultipleFields2LatticeMash< Dim< %(ndim)s > >(const FieldListSet< Dim< %(ndim)s > >& fieldListSet,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                           const FieldList<Dim< %(ndim)s >, int>& mask,
                                           const TableKernel< Dim< %(ndim)s > >& W,
                                           const Dim< %(ndim)s >::Vector& xmin,
                                           const Dim< %(ndim)s >::Vector& xmax,
                                           const vector<int>& nsample,
                                           vector< vector<Dim< %(ndim)s >::Scalar> >& scalarValues,
                                           vector< vector<Dim< %(ndim)s >::Vector> >& vectorValues,
                                           vector< vector<Dim< %(ndim)s >::Tensor> >& tensorValues,
                                           vector< vector<Dim< %(ndim)s >::SymTensor> >& symTensorValues);

}
"""
