//---------------------------------Spheral++----------------------------------//
// VolumeIntegrationFunctions.
// A set of helper methods to numerically evaluate the appropriate kernel
// normalization.
//
// Created by JMO, Mon Jan  6 16:24:19 PST 2003
//----------------------------------------------------------------------------//
#include "VolumeIntegrationFunctions.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include <math.h>
#include <functional>

namespace Spheral {


//------------------------------------------------------------------------------
// Functor to provide the apppriate integrand.
//------------------------------------------------------------------------------
template<typename Dimension>
double
volumeElement(const double r);

template<>
inline
double
volumeElement< Dim<1> >(const double) {
  return 2.0;
}

template<>
inline
double
volumeElement< Dim<2> >(const double r) {
  return 2.0*M_PI*r;
}

template<>
inline
double
volumeElement< Dim<3> >(const double r) {
  return 4.0*M_PI*r*r;
}

template<typename Dimension, typename KernelType>
struct KernelVolumeElement {
  KernelType mW;
  KernelVolumeElement(const KernelType& W): mW(W) {}
  double operator()(const double r) const {
    return volumeElement<Dimension>(r)*mW(r, 1.0);
  }
};    

//------------------------------------------------------------------------------
// Use the trapezoidal rule to evaluate the volume integral of the given kernel.
//------------------------------------------------------------------------------
template<typename Dimension, typename KernelType>
inline
double simpsonsVolumeIntegral(const KernelType& W,
                              const double rMin,
                              const double rMax,
                              const int numBins) {
  return simpsonsIntegration<KernelVolumeElement<Dimension, KernelType>, double, double>
    (KernelVolumeElement<Dimension, KernelType>(W), rMin, rMax, numBins);
}
  
}
