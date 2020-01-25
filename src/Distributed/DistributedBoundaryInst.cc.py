text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/DistributedBoundary.cc"

namespace Spheral {
  template class DistributedBoundary< Dim< %(ndim)s > >;

  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, int>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, unsigned>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeField(Field<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&) const;

  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeFieldVariableSize(Field<Dim< %(ndim)s >, Dim< %(ndim)s >::FacetedVolume>&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeFieldVariableSize(Field<Dim< %(ndim)s >, std::vector<Dim< %(ndim)s >::Scalar> >&) const;
  template void DistributedBoundary<Dim< %(ndim)s > >::beginExchangeFieldVariableSize(Field<Dim< %(ndim)s >, std::vector<Dim< %(ndim)s >::Vector> >&) const;
}
"""
