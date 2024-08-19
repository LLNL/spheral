//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/TableKernel.cc"
#include "Kernel/BSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class TableKernel< Dim<1> >;
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
  template TableKernel<Dim<%1>>::TableKernel(const %(Wname)s<Dim<%1>>&, const unsigned);
""" % {"Wname" : Wname}

text += """
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class TableKernel< Dim<2> >;
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
  template TableKernel<Dim<%2>>::TableKernel(const %(Wname)s<Dim<%2>>&, const unsigned);
""" % {"Wname" : Wname}

text += """
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TableKernel< Dim<3> >;
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
  template TableKernel<Dim<%3>>::TableKernel(const %(Wname)s<Dim<%3>>&, const unsigned);
""" % {"Wname" : Wname}

text += """
#endif
}