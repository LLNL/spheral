//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VolumeIntegrationFunctions.cc"
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
