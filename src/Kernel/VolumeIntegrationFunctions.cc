//---------------------------------Spheral++----------------------------------//
// VolumeIntegrationFunctions.
// A set of helper methods to numerically evaluate the appropriate kernel
// normalization.
//
// Created by JMO, Mon Jan  6 16:24:19 PST 2003
//----------------------------------------------------------------------------//
#include <math.h>
#include <functional>

#include "VolumeIntegrationFunctions.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Geometry/Dimension.hh"
#include "DBC.hh"

namespace Spheral {
namespace KernelSpace {

using namespace std;
//------------------------------------------------------------------------------
// Functor to provide the apppriate integrand.
//------------------------------------------------------------------------------
template<typename Dimension>
double
volumeElement(const double r);

template<>
double
volumeElement< Dim<1> >(const double r) {
  return 2.0;
}

template<>
double
volumeElement< Dim<2> >(const double r) {
  return 2.0*M_PI*r;
}

template<>
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
double simpsonsVolumeIntegral(const KernelType& W,
                              const double rMin,
                              const double rMax,
                              const int numBins) {
  return simpsonsIntegration<KernelVolumeElement<Dimension, KernelType>, double, double>
    (KernelVolumeElement<Dimension, KernelType>(W), rMin, rMax, numBins);
}
  
}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GaussianKernel.hh"
#include "SincKernel.hh"
#include "NSincPolynomialKernel.hh"
#include "NBSplineKernel.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
template
double simpsonsVolumeIntegral<Dim<1>, GaussianKernel< Dim<1> > >
(const GaussianKernel< Dim<1> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<2>, GaussianKernel< Dim<2> > >
(const GaussianKernel< Dim<2> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<3>, GaussianKernel< Dim<3> > >
(const GaussianKernel< Dim<3> >& W, 
 const double rMin, const double rMax, const int numBins);

//------------------------------------------------------------------------------
template
double simpsonsVolumeIntegral<Dim<1>, SincKernel< Dim<1> > >
(const SincKernel< Dim<1> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<2>, SincKernel< Dim<2> > >
(const SincKernel< Dim<2> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<3>, SincKernel< Dim<3> > >
(const SincKernel< Dim<3> >& W, 
 const double rMin, const double rMax, const int numBins);

//------------------------------------------------------------------------------
template
double simpsonsVolumeIntegral<Dim<1>, NSincPolynomialKernel< Dim<1> > >
(const NSincPolynomialKernel< Dim<1> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<2>, NSincPolynomialKernel< Dim<2> > >
(const NSincPolynomialKernel< Dim<2> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<3>, NSincPolynomialKernel< Dim<3> > >
(const NSincPolynomialKernel< Dim<3> >& W, 
 const double rMin, const double rMax, const int numBins);

//------------------------------------------------------------------------------
template
double simpsonsVolumeIntegral<Dim<1>, NBSplineKernel< Dim<1> > >
(const NBSplineKernel< Dim<1> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<2>, NBSplineKernel< Dim<2> > >
(const NBSplineKernel< Dim<2> >& W, 
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral<Dim<3>, NBSplineKernel< Dim<3> > >
(const NBSplineKernel< Dim<3> >& W, 
 const double rMin, const double rMax, const int numBins);

}
}
