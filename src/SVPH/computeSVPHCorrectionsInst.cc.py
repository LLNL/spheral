text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/computeSVPHCorrections.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

template 
void
computeSVPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                       const TableKernel<Dim< %(ndim)s > >& W,
                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& volume,
                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                       Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                       Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                       Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB);

template 
void
computeSVPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                       const TableKernel<Dim< %(ndim)s > >& W,
                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& volume,
                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                       const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB);

}
"""
