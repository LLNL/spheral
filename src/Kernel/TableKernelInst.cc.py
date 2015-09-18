text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TableKernel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace KernelSpace {
    template class TableKernel< Dim< %(ndim)s > >;
"""

for Wname in ("BSplineKernel",
              "W4SplineKernel",
              "GaussianKernel",
              "SuperGaussianKernel",
              "PiGaussianKernel",
              "HatKernel",
              "SincKernel",
              "NSincPolynomialKernel",
              "QuarticSplineKernel",
              "QuinticSplineKernel",
              "NBSplineKernel",
              "WendlandC4Kernel",
              "WendlandC6Kernel"):
    text += """
    template TableKernel< Dim< %%(ndim)s > >::TableKernel(const %(Wname)s< Dim< %%(ndim)s > >&, const int, const double);
""" % {"Wname" : Wname}

text += """
  }
}
"""
