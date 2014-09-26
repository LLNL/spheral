//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeHullVolumes.cc"

namespace Spheral {
  template void computeHullVolumes(const ConnectivityMap<Dim<1> >&, 
                                   const Dim<1>::Scalar kernelExtent,
                                   const FieldList<Dim<1>, Dim<1>::Vector>&,
                                   const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                   FieldList<Dim<1>, Dim<1>::FacetedVolume>&,
                                   FieldList<Dim<1>, Dim<1>::Scalar>&);
  template void computeHullVolumes(const ConnectivityMap<Dim<2> >&, 
                                   const Dim<2>::Scalar kernelExtent,
                                   const FieldList<Dim<2>, Dim<2>::Vector>&,
                                   const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                   FieldList<Dim<2>, Dim<2>::FacetedVolume>&,
                                   FieldList<Dim<2>, Dim<2>::Scalar>&);
  template void computeHullVolumes(const ConnectivityMap<Dim<3> >&, 
                                   const Dim<3>::Scalar kernelExtent,
                                   const FieldList<Dim<3>, Dim<3>::Vector>&,
                                   const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                   FieldList<Dim<3>, Dim<3>::FacetedVolume>&,
                                   FieldList<Dim<3>, Dim<3>::Scalar>&);
}

