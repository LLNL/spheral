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
  namespace KernelSpace {
"""

for Wname in ("TableKernel",
              "GaussianKernel",
              "SincKernel",
              "NSincPolynomialKernel",
              "NBSplineKernel"):
    text += """
    template
    double simpsonsVolumeIntegral< %%(ndim)s, %(Wname)s< %%(ndim)ss > >
    (const %(Wname)s< %%(ndim)s >& W, 
     const double rMin, const double rMax, const int numBins);

"""

text += """
  }
}
"""
