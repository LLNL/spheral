//------------------------------------------------------------------------------
// A collection of intantiations for SphericalSPH 
//------------------------------------------------------------------------------
#include "Kernel/SphericalKernelOslo.hh"
#include "Geometry/Dimension.hh"

#include "computeSPHSumMassDensity.cc"
#include "correctSPHSumMassDensity.cc"
#include "SPHHydroBase.hh"

namespace Spheral {

template void computeSPHSumMassDensity(const ConnectivityMap<Dim<1>>&, 
                                       const SphericalKernelOslo&, 
                                       const bool,
                                       const FieldList<Dim<1>, Dim<1>::Vector>&,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                       FieldList<Dim<1>, Dim<1>::Scalar>&);

template void correctSPHSumMassDensity(const ConnectivityMap<Dim<1>>&, 
                                       const SphericalKernelOslo&, 
                                       const bool,
                                       const FieldList<Dim<1>, Dim<1>::Vector>&,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                       FieldList<Dim<1>, Dim<1>::Scalar>&);

}
