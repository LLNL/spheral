//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/VolumeIntegrationFunctions.cc"
#include "Kernel/TableKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/NBSplineKernel.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template
double simpsonsVolumeIntegral< Dim<1>, TableKernel< Dim<1> > >
(const TableKernel< Dim<1> >& W,
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral< Dim<1>, GaussianKernel< Dim<1> > >
(const GaussianKernel< Dim<1> >& W,
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral< Dim<1>, NBSplineKernel< Dim<1> > >
(const NBSplineKernel< Dim<1> >& W,
 const double rMin, const double rMax, const int numBins);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
double simpsonsVolumeIntegral< Dim<2>, TableKernel< Dim<2> > >
(const TableKernel< Dim<2> >& W,
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral< Dim<2>, GaussianKernel< Dim<2> > >
(const GaussianKernel< Dim<2> >& W,
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral< Dim<2>, NBSplineKernel< Dim<2> > >
(const NBSplineKernel< Dim<2> >& W,
 const double rMin, const double rMax, const int numBins);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
double simpsonsVolumeIntegral< Dim<3>, TableKernel< Dim<3> > >
(const TableKernel< Dim<3> >& W,
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral< Dim<3>, GaussianKernel< Dim<3> > >
(const GaussianKernel< Dim<3> >& W,
 const double rMin, const double rMax, const int numBins);

template
double simpsonsVolumeIntegral< Dim<3>, NBSplineKernel< Dim<3> > >
(const NBSplineKernel< Dim<3> >& W,
 const double rMin, const double rMax, const int numBins);
#endif
}
