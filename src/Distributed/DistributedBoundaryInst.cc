//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class DistributedBoundary< Dim<1> >;
    template class DistributedBoundary< Dim<2> >;
    template class DistributedBoundary< Dim<3> >;

    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, int>&) const;
    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, unsigned>&) const;
    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, Dim<1>::Scalar>&) const;
    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, Dim<1>::Vector>&) const;
    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, Dim<1>::Tensor>&) const;
    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>&) const;
    template void DistributedBoundary<Dim<1> >::beginExchangeField(FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>&) const;

    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, int>&) const;
    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, unsigned>&) const;
    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, Dim<2>::Scalar>&) const;
    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, Dim<2>::Vector>&) const;
    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, Dim<2>::Tensor>&) const;
    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>&) const;
    template void DistributedBoundary<Dim<2> >::beginExchangeField(FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>&) const;

    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, int>&) const;
    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, unsigned>&) const;
    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, Dim<3>::Scalar>&) const;
    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, Dim<3>::Vector>&) const;
    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, Dim<3>::Tensor>&) const;
    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>&) const;
    template void DistributedBoundary<Dim<3> >::beginExchangeField(FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>&) const;
  }
}
