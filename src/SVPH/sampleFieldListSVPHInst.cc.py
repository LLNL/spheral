text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/sampleFieldListSVPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

using std::vector;

// Scalar
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
sampleFieldListSVPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                            const TableKernel< Dim< %(ndim)s > >& W,
                                            const Mesh<Dim< %(ndim)s > >& mesh,
                                            const bool firstOrderConsistent);

// Vector
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
sampleFieldListSVPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                            const TableKernel< Dim< %(ndim)s > >& W,
                                            const Mesh<Dim< %(ndim)s > >& mesh,
                                            const bool firstOrderConsistent);

// Tensor
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
sampleFieldListSVPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                            const TableKernel< Dim< %(ndim)s > >& W,
                                            const Mesh<Dim< %(ndim)s > >& mesh,
                                            const bool firstOrderConsistent);

// SymTensor
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
sampleFieldListSVPH<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                               const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                               const TableKernel< Dim< %(ndim)s > >& W,
                                               const Mesh<Dim< %(ndim)s > >& mesh,
                                               const bool firstOrderConsistent);
}
"""
