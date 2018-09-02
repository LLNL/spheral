text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TableKernel.cc"
#include "Geometry/Dimension.hh"

#include "BSplineKernel.hh"
#include "W4SplineKernel.hh"
#include "GaussianKernel.hh"
#include "SuperGaussianKernel.hh"
#include "PiGaussianKernel.hh"
#include "SincKernel.hh"
#include "NSincPolynomialKernel.hh"
#include "NBSplineKernel.hh"
#include "HatKernel.hh"
#include "QuarticSplineKernel.hh"
#include "QuinticSplineKernel.hh"
#include "WendlandC2Kernel.hh"
#include "WendlandC4Kernel.hh"
#include "WendlandC6Kernel.hh"
#include "ExpInvKernel.hh"

namespace Spheral {
  template class TableKernel< Dim< %(ndim)s > >;
"""

for Wname in ("BSplineKernel",
              "W4SplineKernel",
              "GaussianKernel",
              "SuperGaussianKernel",
              "PiGaussianKernel",
              "SincKernel",
              "NSincPolynomialKernel",
              "NBSplineKernel",
              "HatKernel",
              "QuarticSplineKernel",
              "QuinticSplineKernel",
              "WendlandC2Kernel",
              "WendlandC4Kernel",
              "WendlandC6Kernel",
              "ExpInvKernel"):
    text += """
    template TableKernel< Dim< %%(ndim)s > >::TableKernel(const %(Wname)s< Dim< %%(ndim)s > >&, const int, const double);
//    template void TableKernel< Dim< %%(ndim)s > >::augment(const %(Wname)s< Dim< %%(ndim)s > >&);
""" % {"Wname" : Wname}

text += """
}
"""
