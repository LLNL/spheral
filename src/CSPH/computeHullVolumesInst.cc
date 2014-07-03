//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeHullVolumes.cc"

namespace Spheral {
  template void computeHullVolumes(const ConnectivityMap<Dim<1> >&, 
                                   const FieldList<Dim<1>, Dim<1>::Vector>&,
                                   FieldList<Dim<1>, Dim<1>::Scalar>&);
  template void computeHullVolumes(const ConnectivityMap<Dim<2> >&, 
                                   const FieldList<Dim<2>, Dim<2>::Vector>&,
                                   FieldList<Dim<2>, Dim<2>::Scalar>&);
  template void computeHullVolumes(const ConnectivityMap<Dim<3> >&, 
                                   const FieldList<Dim<3>, Dim<3>::Vector>&,
                                   FieldList<Dim<3>, Dim<3>::Scalar>&);
}

