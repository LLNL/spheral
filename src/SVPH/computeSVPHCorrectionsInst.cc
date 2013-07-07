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
computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dim<1> >& connectivityMap,
                       const KernelSpace::TableKernel<Dim<1> >& W,
                       const FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& volume,
                       const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position,
                       const FieldSpace::FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                       FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& B,
                       FieldSpace::FieldList<Dim<1>, Dim<1>::Tensor>& gradB);

template 
void
computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dim<2> >& connectivityMap,
                       const KernelSpace::TableKernel<Dim<2> >& W,
                       const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& volume,
                       const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& position,
                       const FieldSpace::FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                       FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& B,
                       FieldSpace::FieldList<Dim<2>, Dim<2>::Tensor>& gradB);

template 
void
computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dim<3> >& connectivityMap,
                       const KernelSpace::TableKernel<Dim<3> >& W,
                       const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& volume,
                       const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& position,
                       const FieldSpace::FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                       FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& B,
                       FieldSpace::FieldList<Dim<3>, Dim<3>::Tensor>& gradB);

}
}
