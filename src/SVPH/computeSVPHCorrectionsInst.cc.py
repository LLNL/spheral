text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computeSVPHCorrections.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace SVPHSpace {

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;

template 
void
computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                       const KernelSpace::TableKernel<Dim< %(ndim)s > >& W,
                       const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& volume,
                       const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                       const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                       FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                       FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                       FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB);

template 
void
computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                       const KernelSpace::TableKernel<Dim< %(ndim)s > >& W,
                       const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& volume,
                       const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                       const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                       FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                       FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                       FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB);

}
}
"""
