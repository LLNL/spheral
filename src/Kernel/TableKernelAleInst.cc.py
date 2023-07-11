text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/TableKernel.cc"
#include "Geometry/Dimension.hh"

#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/NBSplineKernel.hh"

namespace Spheral {
  template class TableKernel< Dim< %(ndim)s > >;
"""

for Wname in ("WendlandC2Kernel",
              "WendlandC4Kernel",
              "WendlandC6Kernel",
              "GaussianKernel",
              "PiGaussianKernel",
              "NBSplineKernel"):
    text += """
  template TableKernel<Dim<%%(ndim)s>>::TableKernel(const %(Wname)s<Dim<%%(ndim)s>>&, const unsigned);
""" % {"Wname" : Wname}

text += """
}
"""
