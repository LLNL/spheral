//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "SolidSPH/NodeCoupling.hh"
#include "computeSolidCRKSPHSumMassDensity.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeSolidCRKSPHSumMassDensity(const ConnectivityMap<Dim<1> >&, 
                                                   const TableKernel<Dim<1> >&, 
                                                   const FieldList<Dim<1>, Dim<1>::Vector>&,
                                                   const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                                   const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                                   const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                                   const NodeCoupling&,
                                                   FieldList<Dim<1>, Dim<1>::Scalar>&);
    template void computeSolidCRKSPHSumMassDensity(const ConnectivityMap<Dim<2> >&, 
                                                   const TableKernel<Dim<2> >&, 
                                                   const FieldList<Dim<2>, Dim<2>::Vector>&,
                                                   const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                                   const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                                   const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                                   const NodeCoupling&,
                                                   FieldList<Dim<2>, Dim<2>::Scalar>&);
    template void computeSolidCRKSPHSumMassDensity(const ConnectivityMap<Dim<3> >&, 
                                                   const TableKernel<Dim<3> >&, 
                                                   const FieldList<Dim<3>, Dim<3>::Vector>&,
                                                   const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                                   const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                                   const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                                   const NodeCoupling&,
                                                   FieldList<Dim<3>, Dim<3>::Scalar>&);
  }
}

