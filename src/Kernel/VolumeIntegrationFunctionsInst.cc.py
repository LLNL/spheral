text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VolumeIntegrationFunctions.cc"
#include "TableKernel.hh"
#include "GaussianKernel.hh"
#include "SincKernel.hh"
#include "NSincPolynomialKernel.hh"
#include "NBSplineKernel.hh"

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
