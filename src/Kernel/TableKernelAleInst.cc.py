text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/TableKernel.cc"
#include "Geometry/Dimension.hh"

#include "Kernel/BSplineKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/NBSplineKernel.hh"

namespace Spheral {
  template class TableKernel< Dim< %(ndim)s > >;
"""

for Wname in ("BSplineKernel",
              "QuarticSplineKernel",
              "QuinticSplineKernel",
              "GaussianKernel",
              "PiGaussianKernel",
              "NBSplineKernel"):
    text += """
  template TableKernel<Dim<%%(ndim)s>>::TableKernel(const %(Wname)s<Dim<%%(ndim)s>>&, const unsigned);
""" % {"Wname" : Wname}

text += """
}
"""
