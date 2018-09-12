text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/VolumeIntegrationFunctions.cc"
#include "Kernel/TableKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/NBSplineKernel.hh"

namespace Spheral {
"""

for Wname in ("TableKernel",
              "GaussianKernel",
              "SincKernel",
              "NSincPolynomialKernel",
              "NBSplineKernel"):
    text += """
    template
    double simpsonsVolumeIntegral< Dim< %%(ndim)s >, %(Wname)s< Dim< %%(ndim)s > > >
    (const %(Wname)s< Dim< %%(ndim)s > >& W, 
     const double rMin, const double rMax, const int numBins);

""" % {"Wname" : Wname}

text += """
}
"""
