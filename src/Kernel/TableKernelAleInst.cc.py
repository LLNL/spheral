text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/TableKernel.cc"
#include "Geometry/Dimension.hh"

#include "Kernel/GaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/NBSplineKernel.hh"

namespace Spheral {
  template class TableKernel< Dim< %(ndim)s > >;
"""

for Wname in ("GaussianKernel",
              "PiGaussianKernel",
              "NBSplineKernel"):
    text += """
    template TableKernel< Dim< %%(ndim)s > >::TableKernel(const %(Wname)s< Dim< %%(ndim)s > >&, const int, const double);
//    template void TableKernel< Dim< %%(ndim)s > >::augment(const %(Wname)s< Dim< %%(ndim)s > >&);
""" % {"Wname" : Wname}

text += """
}
"""
