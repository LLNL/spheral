#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 1D objects as generic names.
# This version is specialized for the 1D spherical formulation.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if "1d" in x]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("1d", ""), name))
for name in [x for x in Spheral.__dict__ if (x != "Spherical" and x[:9] == "Spherical")]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("Spherical", ""), name))
del x, name
from Spheral import *
GeometryRegistrar.coords(CoordinateType.Spherical)

# Make our aliases for the kernels point to the 3D instances, since that's what we need for building SphericalKernel
TableKernel = TableKernel3d
BSplineKernel = BSplineKernel3d
NBSplineKernel = NBSplineKernel3d
W4SplineKernel = W4SplineKernel3d
GaussianKernel = GaussianKernel3d
SuperGaussianKernel = SuperGaussianKernel3d
PiGaussianKernel = PiGaussianKernel3d
HatKernel = HatKernel3d
SincKernel = SincKernel3d
NSincPolynomialKernel = NSincPolynomialKernel3d
QuarticSplineKernel = QuarticSplineKernel3d
QuinticSplineKernel = QuinticSplineKernel3d
WendlandC2Kernel = WendlandC2Kernel3d
WendlandC4Kernel = WendlandC4Kernel3d
WendlandC6Kernel = WendlandC6Kernel3d
