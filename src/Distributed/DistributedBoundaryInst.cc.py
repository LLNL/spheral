text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class DistributedBoundary< Dim< %(ndim)s > >;

    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, int>&) const;
    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, unsigned>&) const;
    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&) const;
    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&) const;
    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&) const;
    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&) const;
    template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(FieldSpace::Field<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&) const;
  }
}
"""
