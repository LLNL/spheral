"""
Spheral Kernel module.

Provide the standard SPH/CRK interpolation kernels.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Kernel/BSplineKernel.hh"',
            '"Kernel/W4SplineKernel.hh"',
            '"Kernel/GaussianKernel.hh"',
            '"Kernel/SuperGaussianKernel.hh"',
            '"Kernel/PiGaussianKernel.hh"',
            '"Kernel/HatKernel.hh"',
            '"Kernel/SincKernel.hh"',
            '"Kernel/NSincPolynomialKernel.hh"',
            '"Kernel/NBSplineKernel.hh"',
            '"Kernel/QuarticSplineKernel.hh"',
            '"Kernel/QuinticSplineKernel.hh"',
            '"Kernel/TableKernel.hh"',
            '"Kernel/WendlandC2Kernel.hh"',
            '"Kernel/WendlandC4Kernel.hh"',
            '"Kernel/WendlandC6Kernel.hh"',
            '"Kernel/ExpInvKernel.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Define a preamble function to expose the protected methods of Kernel.
#-------------------------------------------------------------------------------
preamble = """
namespace Spheral {
  template<typename Dimension, typename Descendant>
  class KernelPublicist: public Kernel<Dimension, Descendant> {
  public:
    using Kernel<Dimension, Descendant>::setVolumeNormalization;
    using Kernel<Dimension, Descendant>::setKernelExtent;
    using Kernel<Dimension, Descendant>::setInflectionPoint;
  };
}
"""

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from Kernel import *

for ndim in dims:
    for KT in ("BSpline","W4Spline", "Gaussian", "SuperGaussian", "PiGaussian",
               "Hat", "Sinc", "NSincPolynomial", "NBSpline", "QuarticSpline",
               "QuinticSpline", "Table", "WendlandC2", "WendlandC4", "WendlandC6", "ExpInv"):
        exec('''
%(KT)sKernel%(ndim)id = PYB11TemplateClass(%(KT)sKernel,
                                           template_parameters = ("Dim<%(ndim)i>", "%(KT)sKernel<Dim<%(ndim)i>>"),
                                            cppname = "%(KT)sKernel<Dim<%(ndim)i>>")
''' % {"ndim" : ndim,
       "KT"   : KT})
