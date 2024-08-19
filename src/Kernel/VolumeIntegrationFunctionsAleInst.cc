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
"""

for Wname in ("TableKernel",
              "GaussianKernel",
              "NBSplineKernel"):
    text += """
    template
    double simpsonsVolumeIntegral< Dim< %1 >, %(Wname)s< Dim< %1 > > >
    (const %(Wname)s< Dim< %1 > >& W,
     const double rMin, const double rMax, const int numBins);

""" % {"Wname" : Wname}

text += """
#endif

#if defined(SPHERAL_ENABLE_2D)
"""

for Wname in ("TableKernel",
              "GaussianKernel",
              "NBSplineKernel"):
    text += """
    template
    double simpsonsVolumeIntegral< Dim< %2 >, %(Wname)s< Dim< %2 > > >
    (const %(Wname)s< Dim< %2 > >& W,
     const double rMin, const double rMax, const int numBins);

""" % {"Wname" : Wname}

text += """
#endif

#if defined(SPHERAL_ENABLE_3D)
"""

for Wname in ("TableKernel",
              "GaussianKernel",
              "NBSplineKernel"):
    text += """
    template
    double simpsonsVolumeIntegral< Dim< %3 >, %(Wname)s< Dim< %3 > > >
    (const %(Wname)s< Dim< %3 > >& W,
     const double rMin, const double rMax, const int numBins);

""" % {"Wname" : Wname}

text += """
#endif
}