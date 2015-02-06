//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeHVolumes.cc"

namespace Spheral {
  template void computeHVolumes(const Dim<1>::Scalar kernelExtent,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                FieldList<Dim<1>, Dim<1>::Scalar>&);
  template void computeHVolumes(const Dim<2>::Scalar kernelExtent,
                                const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                FieldList<Dim<2>, Dim<2>::Scalar>&);
  template void computeHVolumes(const Dim<3>::Scalar kernelExtent,
                                const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                FieldList<Dim<3>, Dim<3>::Scalar>&);
}

